import sys
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    get_sun,
    get_body,
    GeocentricTrueEcliptic
)
import datetime

def get_moon_label():
    """
    Returns ONE label describing the Moon's current state, chosen by priority:
      1) 'lunar_eclipse' (a.k.a. "blood moon")
      2) 'blue_moon'     (second full moon in a single calendar month)
      3) 'black_moon'    (second new moon in a single calendar month)
      4) 'harvest_moon'
      5) 'supermoon'
      6) 'micro_moon'
      7) else one of the 8 standard phases:
         'new_moon', 'waxing_crescent', 'first_quarter', 'waxing_gibbous',
         'full_moon', 'waning_gibbous', 'last_quarter', 'waning_crescent'.

    Real astronomy for eclipses, equinox-based harvest moons,
    etc. can be more complex.
    """

    debug_lines = []

    # Location: Amsterdam, Netherlands
    nl_location = EarthLocation(lat=52.3676*u.deg, lon=4.9041*u.deg, height=0*u.m)
    debug_lines.append(f"Using Netherlands (Amsterdam approx) location for topocentric coords.")

    # Current UTC time
    now = Time.now()
    debug_lines.append(f"Current UTC (ISO): {now.isot}")

    # Geocentric positions of Sun and Moon
    moon_pos = get_body('moon', now)
    sun_pos = get_sun(now)

    # 1) PHASE FRACTION = (1 + cos angular separation) / 2
    elongation = moon_pos.separation(sun_pos)
    phase_fraction = (1 + np.cos(elongation.to(u.rad))) / 2  # range 0..1
    debug_lines.append(f"Elongation: {elongation:.3f}, Phase fraction: {phase_fraction:.4f}")

    # 2) DISTANCE from observer to the Moon (km)
    distance_km = moon_pos.distance.to(u.km).value
    debug_lines.append(f"Observer-to-Moon distance (km): {distance_km:.1f}")

    # 3) Approx ecliptic latitude for rough eclipse check, from topocentric coords
    moon_ecliptic = moon_pos.transform_to(GeocentricTrueEcliptic(equinox=now))
    ecliptic_lat_deg = moon_ecliptic.lat.to(u.deg).value
    debug_lines.append(f"Ecliptic lat (deg): {ecliptic_lat_deg:.3f}")

    # Define thresholds
    FULL_MOON_THRESHOLD = 0.95
    NEW_MOON_THRESHOLD  = 0.05
    SUPERMOON_PERIGEE   = 360000
    MICRO_MOON_APOGEE   = 405000
    ECLIPTIC_LAT_THRESH = 1.5

    is_near_full  = (phase_fraction > FULL_MOON_THRESHOLD)
    is_near_new   = (phase_fraction < NEW_MOON_THRESHOLD)
    is_near_eclip = (abs(ecliptic_lat_deg) < ECLIPTIC_LAT_THRESH)

    debug_lines.append(f"is_near_full={is_near_full}, is_near_new={is_near_new}, is_near_eclip={is_near_eclip}")

    # 4) Rough Eclipse Check: near full + near ecliptic
    might_be_eclipse = is_near_full and is_near_eclip

    # 5) Supermoon / Micro Moon if near full
    is_supermoon  = (distance_km < SUPERMOON_PERIGEE and is_near_full)
    is_micro_moon = (distance_km > MICRO_MOON_APOGEE and is_near_full)

    # 6) Harvest Moon: near full in late Sept or early Oct
    d_utc = now.to_datetime()  # Python datetime in UTC
    is_harvest_time = (d_utc.month == 9 or d_utc.month == 10)
    is_harvest_moon = (is_harvest_time and is_near_full)

    # 7) Blue Moon: second full moon in the same calendar month
    is_blue_moon = False
    if is_near_full:
        current_month = d_utc.month
        day_step = 0.5
        t_search = now
        found_another_full = False
        for _ in range(60):  # up to 30 days
            t_search = t_search - day_step*u.day
            dt_utc = t_search.to_datetime()
            if dt_utc.month != current_month:
                break
            mpos = get_body('moon', t_search)
            spos = get_sun(t_search)
            el = mpos.separation(spos)
            frac = (1 + np.cos(el.to(u.rad))) / 2
            if frac > FULL_MOON_THRESHOLD:
                found_another_full = True
                break
        is_blue_moon = found_another_full
    debug_lines.append(f"is_blue_moon={is_blue_moon}")

    # 8) Black Moon: second new moon in the same month
    is_black_moon = False
    if is_near_new:
        current_month = d_utc.month
        day_step = 0.5
        t_search = now
        found_another_new = False
        for _ in range(60):
            t_search = t_search - day_step*u.day
            dt_utc = t_search.to_datetime()
            if dt_utc.month != current_month:
                break
            mpos = get_body('moon', t_search)
            spos = get_sun(t_search)
            el = mpos.separation(spos)
            frac = (1 + np.cos(el.to(u.rad))) / 2
            if frac < NEW_MOON_THRESHOLD:
                found_another_new = True
                break
        is_black_moon = found_another_new
    debug_lines.append(f"is_black_moon={is_black_moon}")

    # 9) Decide final label by priority
    label = None
    if might_be_eclipse:
        label = "lunar_eclipse"
    elif is_blue_moon:
        label = "blue_moon"
    elif is_black_moon:
        label = "black_moon"
    elif is_harvest_moon:
        label = "harvest_moon"
    elif is_supermoon:
        label = "supermoon"
    elif is_micro_moon:
        label = "micro_moon"
    else:
        # Fall back to 8-phase logic
        if is_near_new:
            label = "new_moon"
        elif phase_fraction < 0.1875:
            label = "waxing_crescent"
        elif phase_fraction < 0.3125:
            label = "first_quarter"
        elif phase_fraction < 0.4375:
            label = "waxing_gibbous"
        elif phase_fraction < 0.5625:
            label = "full_moon"
        elif phase_fraction < 0.6875:
            label = "waning_gibbous"
        elif phase_fraction < 0.8125:
            label = "last_quarter"
        elif phase_fraction < 0.9375:
            label = "waning_crescent"
        else:
            label = "new_moon"

    debug_lines.append(f"Final label: {label}")

    # Print debug info to stderr so it doesn't conflict with GH Actions output
    for dbg in debug_lines:
        print(f"[DEBUG] {dbg}", file=sys.stderr)

    return label

if __name__ == "__main__":
    # Output only the final label to stdout
    print(get_moon_label(), end="")
