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

from PIL import Image, ImageFilter

# Shapely for geometry and affine transforms.
from shapely.geometry import Polygon
from matplotlib.patches import PathPatch, Circle
from matplotlib.path import Path

# Astropy for Moon/Sun positions, distances, etc.
from astropy.time import Time
from astropy.coordinates import get_body, get_sun, EarthLocation, AltAz
import astropy.units as u

# Configure logging (human‑readable messages).
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# --- Observer Location (example: Amsterdam, Netherlands) ---
nl_location = EarthLocation(lat=52.3676 * u.deg, lon=4.9041 * u.deg, height=0 * u.m)

# --- Helper Functions ---

def create_circle(center, radius, resolution=512):
    """Return a Shapely Polygon approximating a circle."""
    angles = np.linspace(0, 2 * np.pi, resolution, endpoint=True)
    coords = [(center[0] + radius * np.cos(a),
               center[1] + radius * np.sin(a)) for a in angles]
    return Polygon(coords)

def shapely_poly_to_patch(poly, **kwargs):
    """Convert a Shapely Polygon (or MultiPolygon) into a Matplotlib PathPatch."""
    if poly.is_empty:
        return None
    if poly.geom_type == 'MultiPolygon':
        patches = [shapely_poly_to_patch(subpoly, **kwargs)
                   for subpoly in poly.geoms if not subpoly.is_empty]
        return patches
    coords = np.array(poly.exterior.coords)
    path = Path(coords)
    return PathPatch(path, **kwargs)

def compute_phase_and_position_angle(time_str='2025-02-11T00:00:00', location=None):
    """
    Compute the Moon’s phase data:
      - phase_angle: the Sun–Moon–Earth angle (in degrees). Note that by definition,
        a phase_angle near 0° corresponds to a fully illuminated (full) moon,
        while ~180° means nearly dark (new moon).
      - f: illuminated fraction = (1 + cos(phase_angle))/2.
      - rot_angle: the rotation angle for the terminator so that when the bright‑limb
        is at 90°, the terminator is vertical.
      - elongation: the angular separation (in degrees) between the Moon and the Sun.
      - distance: the geocentric distance to the Moon (in kilometers).

    When a location is provided the Moon’s topocentric coordinates are used;
    otherwise geocentric coordinates are used.
    """
    t = Time(time_str)
    sun = get_sun(t)
    if location is not None:
        moon = get_body('moon', t, location=location)
        altaz_frame = AltAz(obstime=t, location=location)
        moon_altaz = moon.transform_to(altaz_frame)
        sun_altaz = sun.transform_to(altaz_frame)
        pa = np.arctan2(np.sin(sun_altaz.az.radian - moon_altaz.az.radian) *
                        np.cos(sun_altaz.alt.radian),
                        np.sin(sun_altaz.alt.radian) -
                        np.sin(moon_altaz.alt.radian) * np.cos(sun_altaz.az.radian - moon_altaz.az.radian))
        pa_deg = np.degrees(pa) % 360.
        logging.info(f"Topocentric bright‑limb position angle: {pa_deg:.1f}°")
    else:
        moon = get_body('moon', t)
        ra_moon, dec_moon = moon.ra.radian, moon.dec.radian
        ra_sun, dec_sun = sun.ra.radian, sun.dec.radian
        pa = np.arctan2(np.cos(dec_sun) * np.sin(ra_sun - ra_moon),
                        np.sin(dec_sun) * np.cos(dec_moon) -
                        np.cos(dec_sun) * np.sin(dec_moon) * np.cos(ra_sun - ra_moon))
        pa_deg = np.degrees(pa) % 360.
        logging.info(f"Geocentric bright‑limb position angle: {pa_deg:.1f}°")

    # Compute the phase angle using the cosine law.
    d_moon = moon.distance.to(u.km).value
    d_sun = sun.distance.to(u.km).value
    psi = sun.separation(moon).radian  # angular separation (radians)
    r = np.sqrt(d_sun**2 + d_moon**2 - 2 * d_sun * d_moon * np.cos(psi))
    cos_i = (d_moon**2 + r**2 - d_sun**2) / (2 * d_moon * r)
    cos_i = np.clip(cos_i, -1, 1)
    i_angle = np.arccos(cos_i)
    phase_angle = np.degrees(i_angle)
    f = (1 + np.cos(i_angle)) / 2
    rot_angle = -(pa_deg - 90)

    # Compute elongation and distance.
    elongation = sun.separation(moon).degree
    distance = d_moon
    return phase_angle, f, rot_angle, elongation, distance

def solve_d(f, R=1.0, tol=1e-5):
    """
    For two circles (each of radius R) separated by distance d, the area of intersection is:
      A(d) = 2 R² arccos(d/(2R)) - (d/2)*√(4R² - d²)
    We solve for d such that A(d) = π R² · f.
    """
    target = np.pi * R**2 * f
    low, high = 0.0, 2 * R
    while high - low > tol:
        mid = (low + high) / 2
        A_mid = 2 * R**2 * np.arccos(mid/(2*R)) - (mid/2) * np.sqrt(4*R**2 - mid**2)
        if A_mid > target:
            low = mid
        else:
            high = mid
    return (low + high) / 2

def compute_moon_shadow_astropy(time_str, R=1.0, resolution=512, location=None):
    """
    Compute the dark (shadow) region on the Moon’s disk using the two‑circle method.
    This version uses the determined phase stage to choose the auxiliary circle's
    position. For waxing phases the auxiliary circle is centered at (+d, 0) and for
    waning phases at (–d, 0).
    """
    phase_angle_deg, f, rot_angle, elongation, distance = compute_phase_and_position_angle(time_str, location=location)
    white = create_circle((0, 0), R, resolution)

    # For nearly full or nearly new, we use the extremes.
    if f >= 0.999:
        logging.info("Moon is nearly full; using whole disk as bright region.")
        bright = white
        dark = Polygon()  # empty
        return white, dark, phase_angle_deg, f, rot_angle
    if f <= 0.001:
        logging.info("Moon is nearly new; using whole disk as dark region.")
        bright = Polygon()  # empty
        dark = white
        return white, dark, phase_angle_deg, f, rot_angle

    # Solve for the separation d that gives the correct overlapping area.
    d = solve_d(f, R=R)

    # Determine the phase stage (e.g., Waxing Crescent, Waxing Gibbous, etc.)
    stage = determine_moon_stage(phase_angle_deg, f, rot_angle)

    # Choose the auxiliary circle's center based on the stage.
    if stage in ("Waxing Crescent", "Waxing Gibbous"):
        # For waxing, assume the Sun is on the right.
        aux_center = (d, 0)
        # For a crescent (f < 0.5), use intersection to get a thin bright region.
        # For gibbous (f > 0.5), use difference so that only a bite is dark.
        if f < 0.5:
            bright = white.intersection(create_circle(aux_center, R, resolution))
            logging.info("Waxing Crescent: bright region computed via intersection.")
        else:
            bright = white.difference(create_circle(aux_center, R, resolution))
            logging.info("Waxing Gibbous: bright region computed via difference.")
    else:
        # For waning phases, assume the Sun is on the left.
        aux_center = (-d, 0)
        if f < 0.5:
            bright = white.intersection(create_circle(aux_center, R, resolution))
            logging.info("Waning Crescent: bright region computed via intersection.")
        else:
            bright = white.difference(create_circle(aux_center, R, resolution))
            logging.info("Waning Gibbous: bright region computed via difference.")

    dark = white.difference(bright)
    return white, dark, phase_angle_deg, f, rot_angle

def determine_moon_stage(phase_angle, f, rot_angle, tol_f=0.05):
    """
    Return a rough textual stage of the Moon based on the illuminated fraction f
    and the sign of rot_angle (which distinguishes waxing from waning).
      - New Moon when f is very low.
      - Full Moon when f is nearly 1.
      - For intermediate f, we distinguish Crescent vs. Gibbous.
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
    A very simple heuristic to flag a special Moon.
      - For example, if a Full Moon occurs on a day later in the month, mark it as a Blue Moon.
      - (Extend this function with eclipse logic for a Blood Moon if desired.)
    """
    t = Time(time_str)
    day = t.datetime.day
    if stage == "Full Moon" and day >= 20:
        return "Blue Moon"
    return "Normal"

def get_moon_data(time_str, R=1.0, resolution=512, location=None):
    """
    Gather all the Moon data into a dictionary.
    This includes:
      - The white disk and shadow polygons.
      - The phase angle, illuminated fraction, terminator rotation angle.
      - Elongation (Moon–Sun separation in degrees) and Moon distance (km).
      - A textual stage (e.g. 'Waxing Crescent', 'Full Moon', etc.).
      - An approximate moon age (in days) based on the phase.
      - A status (e.g., 'Super Moon' if the Moon is unusually close).
      - Special flags (e.g., Blue Moon, Blood Moon, etc.).
    """
    # Get geometric data for plotting.
    white, dark, phase_angle, f, rot_angle = compute_moon_shadow_astropy(time_str, R, resolution, location)
    # Get additional data (elongation and distance).
    phase_angle2, f2, rot_angle2, elongation, distance = compute_phase_and_position_angle(time_str, location=location)
    # (phase_angle, f, rot_angle) should be consistent with the above.

    stage = determine_moon_stage(phase_angle, f, rot_angle)

    # Compute an approximate Moon age (in days) using a 29.53‑day synodic period.
    synodic_period = 29.53
    if stage == "New Moon":
        moon_age = 0
    elif stage == "Full Moon":
        moon_age = synodic_period / 2
    else:
        # Use the sign of rot_angle to decide waxing (rot_angle < 0) vs. waning.
        if rot_angle < 0:
            # Waxing: phase_angle goes from ~180 (new) down to 0 (full)
            moon_age = (180 - phase_angle) / 360 * synodic_period
        else:
            # Waning: phase_angle goes from 0 (full) up to 180 (new)
            moon_age = (180 + phase_angle) / 360 * synodic_period

    # Determine if the Moon is unusually close or far.
    # (Thresholds are illustrative; adjust as desired.)
    if distance < 370000:
        distance_status = "Super Moon"
    elif distance > 400000:
        distance_status = "Micro Moon"
    else:
        distance_status = "Normal Moon"

    special = determine_special_moon(time_str, stage)

    logging.info(f"Moon data collected:\n"
                 f"  Phase angle: {phase_angle:.1f}°\n"
                 f"  Illuminated fraction: {f:.3f}\n"
                 f"  Rotation angle: {rot_angle:.1f}°\n"
                 f"  Elongation: {elongation:.1f}°\n"
                 f"  Distance: {distance:,.0f} km ({distance_status})\n"
                 f"  Stage: {stage}\n"
                 f"  Moon age: {moon_age:.1f} days\n"
                 f"  Special: {special}")
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
        'white': white,
        'dark': dark,
        'R': R,
        'time_str': time_str
    }

# --- Plot Functions ---

def plot_white_moon(data, filename="output/moon.png"):
    """
    Plot a white Moon disk (white circle with a black border) on a transparent background.
    """
    R = data['R']
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    moon_circle = Circle((0, 0), R, facecolor='white', edgecolor='black', lw=2)
    ax.add_patch(moon_circle)
    ax.set_xlim(-R*1.1, R*1.1)
    ax.set_ylim(-R*1.1, R*1.1)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.savefig(filename, transparent=True, bbox_inches='tight')
    plt.close()
    logging.info(f"White moon image saved to {filename}")

def plot_shadow(data, filename="output/shadow.png"):
    """
    Plot the Moon’s shadow (dark region) as a black shape with 20% opacity.
    """
    R = data['R']
    dark_poly = data['dark']
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    shadow_patch = shapely_poly_to_patch(dark_poly, facecolor='black', edgecolor='none', alpha=0.2)
    if shadow_patch is not None:
        if isinstance(shadow_patch, list):
            for sp in shadow_patch:
                ax.add_patch(sp)
        else:
            ax.add_patch(shadow_patch)
    ax.set_xlim(-R*1.1, R*1.1)
    ax.set_ylim(-R*1.1, R*1.1)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.savefig(filename, transparent=True, bbox_inches='tight')
    plt.close()
    logging.info(f"Shadow image saved to {filename}")

def plot_final(data, filename="output/final.png"):
    """
    Composite the pre‑generated white Moon image and shadow overlay.
    The shadow is rotated by the computed rotation angle and blurred
    (to mimic the soft terminator).
    """
    moon_img = Image.open("output/moon.png").convert("RGBA")
    shadow_img = Image.open("output/shadow.png").convert("RGBA")

    rot_angle = data['rot_angle']
    shadow_img_rot = shadow_img.rotate(rot_angle, resample=Image.BICUBIC, expand=False)

    blur_radius = 5  # Adjust as desired.
    shadow_img_rot = shadow_img_rot.filter(ImageFilter.GaussianBlur(radius=blur_radius))

    final_img = Image.alpha_composite(moon_img, shadow_img_rot)
    final_img.save(filename)
    logging.info(f"Final composite image saved to {filename}")

# --- Main Execution ---
if __name__ == '__main__':
    logging.info("Starting moon phase image generation with extended data")
    os.makedirs("output", exist_ok=True)
    # Set the observation time (UTC).
    now = Time.now()
    time_str = now.to_datetime()
    time_str = '2025-03-02T00:00:00'
    logging.info(f"Observation time: {time_str}")

    R = 1  # Moon radius in arbitrary units.

    # Gather all moon data (using our observer location for phase details).
    moon_data = get_moon_data(time_str, R, resolution=512, location=nl_location)

    # Generate images using the collected data.
    plot_white_moon(moon_data, filename="output/moon.png")
    plot_shadow(moon_data, filename="output/shadow.png")
    plot_final(moon_data, filename="output/final.png")

    logging.info("Moon phase image generation complete")
