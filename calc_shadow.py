#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
import logging
import warnings
from datetime import datetime

from astropy.utils import iers
iers.conf.auto_download = False

# Suppress non‑rotation transformation warnings.
from astropy.coordinates.baseframe import NonRotationTransformationWarning
warnings.filterwarnings("ignore", category=NonRotationTransformationWarning)

from astropy.time import Time
from astropy.coordinates import get_body, get_sun, EarthLocation, AltAz
import astropy.units as u

from PIL import Image, ImageFilter

# Configure logging (human‑readable messages).
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# --- Observer Location (example: Amsterdam, Netherlands) ---
nl_location = EarthLocation(lat=52.3676*u.deg, lon=4.9041*u.deg, height=0*u.m)

# --- Ephemeris and Moon Data Computations ---

def compute_phase_and_position_angle(time_str='2025-02-11T00:00:00', location=None):
    """
    Compute the Moon’s phase data:
      - phase_angle: Sun–Moon–Earth angle (degrees). By definition, near 0° means nearly full,
        near 180° means nearly new.
      - f: illuminated fraction = (1 + cos(phase_angle))/2.
      - rot_angle: rotation angle for the terminator so that when the bright‑limb is at 90°
        the terminator is vertical.
      - elongation: angular separation (degrees) between the Moon and the Sun.
      - distance: geocentric distance to the Moon (km).

    When a location is provided, topocentric coordinates are used.
    """
    t = Time(time_str)
    sun = get_sun(t)
    if location is not None:
        moon = get_body('moon', t, location=location)
        altaz_frame = AltAz(obstime=t, location=location)
        moon_altaz = moon.transform_to(altaz_frame)
        sun_altaz  = sun.transform_to(altaz_frame)
        pa = np.arctan2(np.sin(sun_altaz.az.radian - moon_altaz.az.radian) *
                        np.cos(sun_altaz.alt.radian),
                        np.sin(sun_altaz.alt.radian) -
                        np.sin(moon_altaz.alt.radian) *
                        np.cos(sun_altaz.az.radian - moon_altaz.az.radian))
        pa_deg = np.degrees(pa) % 360.
        logging.info(f"Topocentric bright‑limb position angle: {pa_deg:.1f}°")
    else:
        moon = get_body('moon', t)
        ra_moon, dec_moon = moon.ra.radian, moon.dec.radian
        ra_sun, dec_sun   = sun.ra.radian, sun.dec.radian
        pa = np.arctan2(np.cos(dec_sun)*np.sin(ra_sun - ra_moon),
                        np.sin(dec_sun)*np.cos(dec_moon) -
                        np.cos(dec_sun)*np.sin(dec_moon)*np.cos(ra_sun - ra_moon))
        pa_deg = np.degrees(pa) % 360.
        logging.info(f"Geocentric bright‑limb position angle: {pa_deg:.1f}°")
    # Compute phase angle using the cosine law.
    d_moon = moon.distance.to(u.km).value
    d_sun  = sun.distance.to(u.km).value
    psi    = sun.separation(moon).radian
    r      = np.sqrt(d_sun**2 + d_moon**2 - 2*d_sun*d_moon*np.cos(psi))
    cos_i  = (d_moon**2 + r**2 - d_sun**2) / (2*d_moon*r)
    cos_i  = np.clip(cos_i, -1, 1)
    i_angle = np.arccos(cos_i)
    phase_angle = np.degrees(i_angle)
    f = (1 + np.cos(i_angle)) / 2
    rot_angle = -(pa_deg - 90)
    elongation = sun.separation(moon).degree
    distance = d_moon
    logging.info(f"Phase angle: {phase_angle:.1f}°, f: {f:.3f}, rot_angle: {rot_angle:.1f}°, " +
                 f"elongation: {elongation:.1f}°, distance: {distance:,.0f} km")
    return phase_angle, f, rot_angle, elongation, distance

def determine_moon_stage(phase_angle, f, rot_angle, tol_f=0.05):
    """
    Return a rough textual stage of the Moon based on illuminated fraction f and
    the sign of rot_angle (which distinguishes waxing from waning).
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
    A simple heuristic to flag a special Moon (e.g., Blue Moon).
    """
    t = Time(time_str)
    day = t.datetime.day
    if stage == "Full Moon" and day >= 20:
        return "Blue Moon"
    return "Normal"

def get_moon_data(time_str, R=1.0, resolution=512, location=None):
    """
    Compute and gather all available Moon data into a dictionary.
    This includes:
      - phase angle, illuminated fraction, terminator rotation angle.
      - elongation (Moon–Sun separation) and Moon distance.
      - A textual stage (e.g., "Waxing Crescent", "Full Moon", etc.).
      - An approximate Moon age (days, using a 29.53‑day synodic period).
      - A distance status ("Super Moon", "Micro Moon", "Normal Moon").
      - A special flag (e.g., "Blue Moon").
    """
    phase_angle, f, rot_angle, elongation, distance = compute_phase_and_position_angle(time_str, location=location)
    stage = determine_moon_stage(phase_angle, f, rot_angle)
    synodic_period = 29.53  # days
    if stage == "New Moon":
        moon_age = 0.0
    elif stage == "Full Moon":
        moon_age = synodic_period / 2
    else:
        if rot_angle < 0:  # waxing
            moon_age = (180 - phase_angle) / 360 * synodic_period
        else:              # waning
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

# --- Professional Spherical Terminator Computation (Pixel-by-Pixel) ---

def get_sun_vector(time_str, location=None):
    """
    Compute the unit vector (in Moon‑centered Cartesian coordinates) pointing from the Moon to the Sun.
    We assume the observer's view is along +z so that the visible hemisphere is defined by z ≥ 0.
    """
    t = Time(time_str)
    if location is not None:
        moon = get_body('moon', t, location=location)
    else:
        moon = get_body('moon', t)
    sun = get_sun(t)
    # Get Cartesian coordinates (in km)
    moon_cart = moon.cartesian.xyz.value.flatten()
    sun_cart  = sun.cartesian.xyz.value.flatten()
    S_vec = sun_cart - moon_cart
    S = S_vec / np.linalg.norm(S_vec)
    logging.info(f"Sun vector (Moon-to-Sun): [{S[0]:.3f}, {S[1]:.3f}, {S[2]:.3f}]")
    return S

def plot_shadow(time_str, R=1.0, resolution=512, filename="output/shadow.png", location=None):
    """
    Create an image of the Moon's shadow overlay using a rigorous, spherical, pixel-by-pixel method.

    A grid of (x,y) points covering the Moon's disk (x²+y² ≤ R²) is generated. For each point,
    z = sqrt(R² - x² - y²) is computed (for the visible hemisphere). The illumination is computed as
          illumination = (x,y,z) · S,
    where S is the unit vector from the Moon to the Sun.

    Pixels with illumination > 0 (lit) are transparent; pixels with illumination ≤ 0 (in shadow)
    are rendered as black with 20% opacity.
    """
    S = get_sun_vector(time_str, location=location)
    xs = np.linspace(-R, R, resolution)
    ys = np.linspace(-R, R, resolution)
    X, Y = np.meshgrid(xs, ys)
    R2 = X**2 + Y**2
    mask = R2 <= R**2
    Z = np.zeros_like(X)
    # For points inside the disk, compute z (the positive square root)
    Z[mask] = np.sqrt(R**2 - R2[mask])
    # Compute illumination (dot product with S)
    illumination = np.full_like(X, -1e10, dtype=float)
    illumination[mask] = X[mask]*S[0] + Y[mask]*S[1] + Z[mask]*S[2]
    # Create an alpha channel:
    #   lit pixels (illumination > 0) get alpha 0 (transparent),
    #   dark pixels get alpha 51 (~20% opacity of 255).
    alpha = np.zeros_like(X, dtype=float)
    alpha[mask] = np.where(illumination[mask] > 0, 0, 51)
    # Create an RGBA image array (black color with computed alpha)
    img_array = np.zeros((resolution, resolution, 4), dtype=np.uint8)
    img_array[..., 3] = alpha.astype(np.uint8)
    shadow_img = Image.fromarray(img_array, mode="RGBA")
    shadow_img.save(filename)
    logging.info(f"Shadow image saved to {filename}")

# --- White Moon Rendering ---
def plot_white_moon(R=1.0, resolution=512, filename="output/moon.png"):
    """
    Create an image of the white Moon (a white circle with a black border)
    on a transparent background, forced to resolution x resolution pixels.

    The axis is set so that x and y run exactly from -R to R.
    """
    # Create a figure with the desired pixel resolution.
    fig, ax = plt.subplots(figsize=(resolution/100, resolution/100), dpi=100)
    from matplotlib.patches import Circle
    # Draw the Moon circle (centered at (0,0) with radius R)
    moon_circle = Circle((0, 0), R, facecolor='white', edgecolor='black', lw=2)
    ax.add_patch(moon_circle)
    # Set the coordinate system exactly to [-R, R] in both directions.
    ax.set_xlim(-R, R)
    ax.set_ylim(-R, R)
    ax.set_aspect('equal')
    ax.axis('off')
    # Remove any margins.
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    # Save the figure without any additional padding.
    plt.savefig(filename, transparent=True, bbox_inches=None, pad_inches=0)
    plt.close()


# --- Final Composite ---
def plot_final(filename="output/final.png"):
    """
    Composite the white Moon image and the shadow overlay.
    A slight Gaussian blur is applied to the shadow to simulate a soft terminator.
    """
    moon_img = Image.open("output/moon.png").convert("RGBA")
    shadow_img = Image.open("output/shadow.png").convert("RGBA")
    blur_radius = 2  # Adjust as desired.
    shadow_img = shadow_img.filter(ImageFilter.GaussianBlur(radius=blur_radius))
    final_img = Image.alpha_composite(moon_img, shadow_img)
    final_img.save(filename)
    logging.info(f"Final composite image saved to {filename}")

# --- Main Execution ---
if __name__ == '__main__':
    logging.info("Starting professional moon phase image generation with extended data")
    os.makedirs("output", exist_ok=True)
    # Set the observation time (UTC); adjust as needed.
    time_str = '2025-02-11T00:00:00'
    logging.info(f"Observation time: {time_str}")
    R = 1  # Moon radius in arbitrary units.
    resolution = 512  # Use a fixed resolution for both images.

    # Get all the moon data.
    moon_data = get_moon_data(time_str, R, resolution=resolution, location=nl_location)

    # Generate images.
    plot_white_moon(R, resolution=resolution, filename="output/moon.png")
    plot_shadow(time_str, R, resolution=resolution, filename="output/shadow.png", location=nl_location)
    plot_final(filename="output/final.png")

    logging.info("Moon phase image generation complete")
    # Optionally, log all moon data.
    for key, value in moon_data.items():
        logging.info(f"{key}: {value}")
