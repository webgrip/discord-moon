#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import logging
import warnings
from datetime import datetime

import spiceypy as spice
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u

from PIL import Image, ImageFilter

# Configure logging to show both INFO and DEBUG messages.
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# --- SPICE Kernel Loading ---
# Load the SPICE kernels using your meta-kernel.
# (Make sure the file "moon.tm" is in the current working directory.)
logging.info("Loading SPICE kernels from moon.tm")
spice.furnsh("moon.tm")
# (When done, you can later call spice.kclear())

# --- Observer Location ---
# Example: Amsterdam, Netherlands.
nl_location = EarthLocation(lat=52.3676 * u.deg, lon=4.9041 * u.deg, height=0 * u.m)

def transform_sun_vector_to_local(S, moon_pos):
    """
    Transform the Sun vector S (in J2000) into the Moon's local coordinate system
    where the visible hemisphere is defined by z >= 0. The transformation rotates
    the sub-observer vector (from Moon to Earth) to [0, 0, 1].
    """
    # Compute the sub-observer vector (from Moon to Earth, normalized)
    sub_obs = -np.array(moon_pos) / np.linalg.norm(moon_pos)

    # Define the desired direction (local zenith) as [0, 0, 1]
    target = np.array([0.0, 0.0, 1.0])

    # Compute the rotation axis (cross product) and angle between sub_obs and target.
    axis = np.cross(sub_obs, target)
    axis_norm = np.linalg.norm(axis)

    # If the sub_obs is already aligned (or anti-aligned) with target, no rotation is needed.
    if axis_norm < 1e-6:
        return S  # They are almost colinear.

    axis = axis / axis_norm
    angle = np.arccos(np.clip(np.dot(sub_obs, target), -1.0, 1.0))

    # Rodrigues' rotation formula:
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle)*K + (1 - np.cos(angle))*(K @ K)

    # Transform the Sun vector S into the local frame.
    S_local = R.dot(S)
    return S_local

# --- Functions for Ephemeris & Moon Data Computation using SPICE ---

def compute_phase_and_position_angle(time_str, location):
    """
    Using SPICE, compute:
      - phase_angle: the Sun–Moon–Earth angle in degrees.
      - f: illuminated fraction = (1 + cos(phase_angle))/2.
      - rot_angle: an approximate rotation angle for the terminator (degrees).
      - elongation: angular separation (degrees) between the Moon and Sun as seen from Earth.
      - distance: geocentric distance to the Moon in km.
    """
    # Convert UTC to Ephemeris Time (ET)
    try:
        et = spice.utc2et(time_str)
    except spice.utils.exceptions.SpiceyError as e:
        logging.error("Error converting UTC to ET. Have you loaded the leapseconds kernel?")
        raise

    # Get the Moon state relative to Earth in the J2000 frame.
    # (spkpos returns the position (km) and light time.)
    moon_pos, lt1 = spice.spkpos("MOON", et, "J2000", "NONE", "EARTH")
    # Get the Sun state relative to the Moon.
    sun_pos, lt2 = spice.spkpos("SUN", et, "J2000", "NONE", "MOON")

    # Compute the phase angle: the angle between (Moon->Earth) and (Moon->Sun)
    moon_to_earth = -np.array(moon_pos)  # from Moon to Earth
    moon_to_sun   = np.array(sun_pos)
    # Normalize vectors
    v_me = moon_to_earth / np.linalg.norm(moon_to_earth)
    v_ms = moon_to_sun / np.linalg.norm(moon_to_sun)
    cos_phase = np.dot(v_me, v_ms)
    phase_angle = np.degrees(np.arccos(np.clip(cos_phase, -1, 1)))
    # Illuminated fraction.
    f = (1 + np.cos(np.radians(phase_angle)))/2
    # An approximate terminator rotation angle (not used for pixel shading but for logging)
    # Here we simply set rot_angle = (90 - phase_angle/2) for demonstration.
    rot_angle = 90 - phase_angle/2

    # Compute elongation: angular separation between Sun and Moon as seen from Earth.
    # Use the dot product formula. (Moon position relative to Earth is moon_pos.)
    # Get Sun state relative to Earth.
    sun_earth, lt3 = spice.spkpos("SUN", et, "J2000", "NONE", "EARTH")
    # Normalize
    v_moon = np.array(moon_pos) / np.linalg.norm(moon_pos)
    v_sun  = np.array(sun_earth) / np.linalg.norm(sun_earth)
    cos_elon = np.dot(v_moon, v_sun)
    elongation = np.degrees(np.arccos(np.clip(cos_elon, -1, 1)))

    # Distance from Earth to Moon (in km)
    distance = np.linalg.norm(moon_pos)

    logging.info(f"Phase angle: {phase_angle:.1f}°, f: {f:.3f}, rot_angle: {rot_angle:.1f}°, " +
                 f"elongation: {elongation:.1f}°, distance: {distance:,.0f} km")
    return phase_angle, f, rot_angle, elongation, distance

def determine_moon_stage(phase_angle, f, rot_angle, tol_f=0.05):
    """
    Return a rough textual stage of the Moon.
    """
    if f <= tol_f:
        return "New Moon"
    if f >= 1 - tol_f:
        return "Full Moon"
    if f < 0.5:
        return "Waxing Crescent" if rot_angle < 0 else "Waning Crescent"
    if f > 0.5:
        return "Waxing Gibbous" if rot_angle < 0 else "Waning Gibbous"
    return "Quarter Moon"

def determine_special_moon(time_str, stage):
    """
    A simple heuristic to flag a special Moon.
    """
    t = Time(time_str)
    day = t.datetime.day
    if stage == "Full Moon" and day >= 20:
        return "Blue Moon"
    return "Normal"

def get_moon_data(time_str, R=1.0, resolution=512, location=None):
    """
    Compute and return a dictionary of Moon data.
    """
    phase_angle, f, rot_angle, elongation, distance = compute_phase_and_position_angle(time_str, location=location)
    stage = determine_moon_stage(phase_angle, f, rot_angle)
    synodic_period = 29.53  # days
    if stage == "New Moon":
        moon_age = 0.0
    elif stage == "Full Moon":
        moon_age = synodic_period / 2
    else:
        if rot_angle < 0:
            moon_age = (180 - phase_angle) / 360 * synodic_period
        else:
            moon_age = (180 + phase_angle) / 360 * synodic_period
    if distance < 370000:
        distance_status = "Super Moon"
    elif distance > 400000:
        distance_status = "Micro Moon"
    else:
        distance_status = "Normal Moon"
    special = determine_special_moon(time_str, stage)
    logging.info("Moon data collected:")
    logging.info(f"  Phase angle: {phase_angle:.1f}°")
    logging.info(f"  Illuminated fraction: {f:.3f}")
    logging.info(f"  Rotation angle: {rot_angle:.1f}°")
    logging.info(f"  Elongation: {elongation:.1f}°")
    logging.info(f"  Distance: {distance:,.0f} km ({distance_status})")
    logging.info(f"  Stage: {stage}")
    logging.info(f"  Moon age: {moon_age:.1f} days")
    logging.info(f"  Special: {special}")
    return {
        'phase_angle': phase_angle,
        'f': f,
        'rot_angle': rot_angle,
        'elongation': elongation,
        'distance': distance,
        'stage': stage,
        'moon_age': moon_age,
        'distance_status': distance_status,
        'special': special,
        'time_str': time_str,
        'location': location
    }

# --- Professional Spherical Terminator Computation using SPICE ---
def get_sun_vector_spice(time_str, location):
    """
    Compute the unit vector (in Moon-centered coordinates) from the Moon to the Sun using SPICE.
    Then, transform this vector into the Moon's local frame (with z pointing toward the observer).
    """
    if location is None:
        raise ValueError("A valid observer location is required for SPICE-based calculations.")
    et = spice.utc2et(time_str)
    # Get the Sun position relative to the Moon in the J2000 frame.
    sun_pos, lt = spice.spkpos("SUN", et, "J2000", "NONE", "MOON")
    sun_vec = np.array(sun_pos)
    S = sun_vec / np.linalg.norm(sun_vec)
    logging.info(f"SPICE Moon-to-Sun vector: {S}")

    # Get the Moon position relative to Earth so we can compute the sub-observer vector.
    moon_pos, lt = spice.spkpos("MOON", et, "J2000", "NONE", "EARTH")

    # Transform S into the Moon's local frame so that the sub-observer direction becomes [0, 0, 1].
    S_local = transform_sun_vector_to_local(S, moon_pos)

    # For a nearly full Moon, force S_local to be [0, 0, 1].
    phase_angle, f, rot_angle, elongation, distance = compute_phase_and_position_angle(time_str, location=location)
#     if f > 0.99:
#         S_final = np.array([0.0, 0.0, 1.0])
#     else:
    S_final = S_local
    logging.info(f"Final SPICE Sun vector for illumination (local frame): {S_final}")
    return S_final


def plot_shadow(moon_data, time_str, R=1.0, resolution=512, filename="output/shadow.png", location=None):
    """
    Create an image of the Moon's shadow overlay using a rigorous, spherical, pixel-by-pixel method.

    A grid of (x, y) points covering the Moon's disk (x²+y² ≤ R²) is generated. For each pixel,
    we compute z = sqrt(R² - x² - y²) (the visible hemisphere). Then the illumination is computed as:

         illumination = (x, y, z) · S,

    where S is the unit Sun vector computed by get_sun_vector_spice.

    Pixels with illumination > 0 (lit) are rendered fully transparent; pixels with illumination ≤ 0 (in shadow)
    are rendered as black with ~20% opacity.
    """
    S = get_sun_vector_spice(time_str, location)

    # Create grid over [-R, R] in x and y.
    xs = np.linspace(-R, R, resolution)
    ys = np.linspace(-R, R, resolution)
    X, Y = np.meshgrid(xs, ys)
    R2 = X**2 + Y**2
    mask = R2 <= R**2
    logging.debug(f"Grid shape: {X.shape}")
    logging.debug(f"Number of pixels inside disk: {np.sum(mask)}")

    # Compute z for pixels inside the disk.
    Z = np.zeros_like(X)
    Z[mask] = np.sqrt(R**2 - R2[mask])

    # Compute illumination: dot product of (x,y,z) with S.
    illumination = np.full_like(X, -1e10, dtype=float)
    illumination[mask] = X[mask]*S[0] + Y[mask]*S[1] + Z[mask]*S[2]

    logging.debug(f"Illumination stats (inside disk): min = {np.min(illumination[mask]):.3f}, "
                  f"max = {np.max(illumination[mask]):.3f}, mean = {np.mean(illumination[mask]):.3f}, "
                  f"median = {np.median(illumination[mask]):.3f}")
    lit_pixels = np.sum(illumination[mask] > 0)
    shadow_pixels = np.sum(illumination[mask] <= 0)
    logging.debug(f"Number of lit pixels: {lit_pixels}, shadow pixels: {shadow_pixels}")
    percent_lit = lit_pixels / (lit_pixels + shadow_pixels) * 100
    logging.info(f"Percentage of lit pixels: {percent_lit:.1f}%")

    if (percent_lit != moon_data['f']*100):
        logging.warning(f"Percentage of lit pixels does not match f: {percent_lit:.1f}% vs {moon_data['f']*100:.1f}%")

    # Create an alpha channel: pixels with illumination > 0 (lit) have alpha 0;
    # pixels with illumination <= 0 (shadow) get alpha 51 (~20% opacity).
    alpha = np.zeros_like(X, dtype=float)
    alpha[mask] = np.where(illumination[mask] > 0, 0, 51)
    logging.debug(f"Alpha channel stats (inside disk): min = {np.min(alpha[mask]):.1f}, "
                  f"max = {np.max(alpha[mask]):.1f}, mean = {np.mean(alpha[mask]):.1f}")

    # Build the RGBA image array (black with the computed alpha).
    img_array = np.zeros((resolution, resolution, 4), dtype=np.uint8)
    img_array[..., 3] = alpha.astype(np.uint8)
    shadow_img = Image.fromarray(img_array, mode="RGBA")
    shadow_img.save(filename)
    logging.info(f"Shadow image saved to {filename}")

def plot_white_moon(R=1.0, resolution=512, filename="output/moon.png"):
    """
    Create an image of the white Moon (a white circle with a black border) on a transparent background,
    forced to resolution x resolution pixels. The coordinate system is set so that x and y run exactly from -R to R.
    """
    fig, ax = plt.subplots(figsize=(resolution/100, resolution/100), dpi=100)
    from matplotlib.patches import Circle
    moon_circle = Circle((0, 0), R, facecolor='white', edgecolor='black', lw=2)
    ax.add_patch(moon_circle)
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.savefig(filename, transparent=True, bbox_inches=None, pad_inches=0)
    plt.close()
    logging.info(f"White moon image saved to {filename}")

def plot_final(filename="output/final.png"):
    """
    Composite the white Moon image and the shadow overlay.
    A slight Gaussian blur is applied to the shadow to simulate a soft terminator.
    """
    moon_img = Image.open("output/moon.png").convert("RGBA")
    shadow_img = Image.open("output/shadow.png").convert("RGBA")
    blur_radius = 2  # Adjust as needed.
    shadow_img = shadow_img.filter(ImageFilter.GaussianBlur(radius=blur_radius))
    final_img = Image.alpha_composite(moon_img, shadow_img)
    final_img.save(filename)
    logging.info(f"Final composite image saved to {filename}")

def main():
    logging.info("Starting professional moon phase image generation with SPICE/SpiceyPy")
    os.makedirs("output", exist_ok=True)
    # Set the observation time (UTC). Adjust as needed.
    time_str = '2025-02-12T00:00:00'
    logging.info(f"Observation time: {time_str}")
    R = 1  # Moon radius in arbitrary units.
    resolution = 512  # Use fixed resolution for images.

    # Compute Moon data.
    moon_data = get_moon_data(time_str, R, resolution=resolution, location=nl_location)

    # Generate the white Moon image.
    plot_white_moon(R, resolution=resolution, filename="output/moon.png")
    # Generate the shadow overlay using SPICE-derived Sun vector.
    plot_shadow(moon_data, time_str, R, resolution=resolution, filename="output/shadow.png", location=nl_location)
    # Composite the images.
    plot_final(filename="output/final.png")

    logging.info("Moon phase image generation complete")
    for key, value in moon_data.items():
        logging.info(f"{key}: {value}")

    # When done, clear all loaded SPICE kernels.
    spice.kclear()

if __name__ == '__main__':
    main()
