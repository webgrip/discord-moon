import sys
import math
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import (EarthLocation, get_sun, get_body,
                                 GeocentricTrueEcliptic)
import datetime

def get_moon_label():
    """
    Returns one label describing the Moon's current state at an observer
    in Amsterdam, using:
      - Real geometry for special checks (eclipse, supermoon, etc.)
      - 29.53-day cycle fraction for naming the 8 standard phases
        (new_moon, waxing_crescent, first_quarter, etc.).

    PRIORITY:
      1) 'lunar_eclipse'
      2) 'blue_moon'
      3) 'black_moon'
      4) 'harvest_moon'
      5) 'supermoon'
      6) 'micro_moon'
      Else => one of the 8 time-based phases.
    """

    debug_lines = []

    # 1) Setup: location & current time
    nl_location = EarthLocation(lat=52.3676*u.deg, lon=4.9041*u.deg, height=0*u.m)
    now = Time.now()
    d_utc = now.to_datetime()

    debug_lines.append(f"Using Netherlands (Amsterdam) location.")
    debug_lines.append(f"Current UTC (ISO): {now.isot}")

    # 2) Real geometry for special checks:
    #    (geocentric positions for distance, ecliptic lat, etc.)
    moon_pos = get_body('moon', now)  # Geocentric coordinates of Moon
    sun_pos  = get_sun(now)           # Geocentric coordinates of Sun

    # 2a) Illumination fraction (for logs) => (1 - cos θ)/2
    elongation = moon_pos.separation(sun_pos)
    illum_fraction = (1 - np.cos(elongation.to(u.rad))) / 2
    debug_lines.append(f"Elongation = {elongation:.3f}, IllumFraction ~ {illum_fraction:.4f}")

    # 2b) Distance to Moon (for super/micro checks)
    distance_km = moon_pos.distance.to(u.km).value
    debug_lines.append(f"Geocentric Moon distance (km) = {distance_km:.1f}")

    # 2c) Ecliptic lat (rough eclipse check)
    moon_ecliptic = moon_pos.transform_to(GeocentricTrueEcliptic(equinox=now))
    ecliptic_lat_deg = moon_ecliptic.lat.to(u.deg).value
    debug_lines.append(f"Ecliptic lat (deg) = {ecliptic_lat_deg:.3f}")

    # 3) Basic thresholds
    FULL_MOON_THRESHOLD = 0.995
    NEW_MOON_THRESHOLD  = 0.05
    SUPERMOON_PERIGEE   = 360000  # km
    MICRO_MOON_APOGEE   = 405000  # km
    ECLIPTIC_LAT_THRESH = 1.5

    is_near_full  = (illum_fraction > FULL_MOON_THRESHOLD)
    is_near_new   = (illum_fraction < NEW_MOON_THRESHOLD)
    is_near_eclip = (abs(ecliptic_lat_deg) < ECLIPTIC_LAT_THRESH)

    debug_lines.append(f"is_near_full={is_near_full}, is_near_new={is_near_new}, is_near_eclip={is_near_eclip}")

    # 3a) Eclipse check: near full + near eclip
    might_be_eclipse = (is_near_full and is_near_eclip)

    # 3b) Supermoon / Micro moon if near full
    is_supermoon  = (distance_km < SUPERMOON_PERIGEE and is_near_full)
    is_micro_moon = (distance_km > MICRO_MOON_APOGEE  and is_near_full)

    # 3c) Harvest moon (rough) => near full in late Sept or early Oct
    is_harvest_time = (d_utc.month == 9 or d_utc.month == 10)
    is_harvest_moon = (is_harvest_time and is_near_full)

    # 4) Blue Moon: second near-full in same month
    is_blue_moon = False
    if is_near_full:
        current_month = d_utc.month
        day_step = 0.5
        t_search = now
        found_another_full = False
        for _ in range(60):  # up to 30 days
            t_search = t_search - (day_step*u.day)
            dt_utc_2 = t_search.to_datetime()
            if dt_utc_2.month != current_month:
                break
            # check if near-full there
            mpos2 = get_body('moon', t_search)
            spos2 = get_sun(t_search)
            el2 = mpos2.separation(spos2)
            frac2 = (1 - np.cos(el2.to(u.rad))) / 2
            if frac2 > FULL_MOON_THRESHOLD:
                found_another_full = True
                break
        is_blue_moon = found_another_full

    debug_lines.append(f"is_blue_moon={is_blue_moon}")

    # 5) Black Moon: second near-new in same month
    is_black_moon = False
    if is_near_new:
        current_month = d_utc.month
        day_step = 0.5
        t_search = now
        found_another_new = False
        for _ in range(60):
            t_search = t_search - (day_step*u.day)
            dt_utc_2 = t_search.to_datetime()
            if dt_utc_2.month != current_month:
                break
            mpos2 = get_body('moon', t_search)
            spos2 = get_sun(t_search)
            el2 = mpos2.separation(spos2)
            frac2 = (1 - np.cos(el2.to(u.rad))) / 2
            if frac2 < NEW_MOON_THRESHOLD:
                found_another_new = True
                break
        is_black_moon = found_another_new

    debug_lines.append(f"is_black_moon={is_black_moon}")

    # 6) TIME-BASED FRACTION of the LUNAR CYCLE (for final 8-phase fallback)
    #    We'll pick a reference new moon near 2000-01-06 18:14 UTC (Julian ~2451550.1)
    #    Then compute how far we are into the 29.53-day cycle.
    epoch_jan_2000 = 2451550.1
    synodic_month = 29.53058867
    days_since_ref = now.jd - epoch_jan_2000
    cycle_fraction = (days_since_ref % synodic_month) / synodic_month
    debug_lines.append(f"Cycle fraction = {cycle_fraction:.4f} (time-based)")

    # --> Compute the Moon's age in days (time since last new moon)
    moon_age_days = cycle_fraction * synodic_month
    debug_lines.append(f"Moon Age ~ {moon_age_days:.2f} days since last new moon")

    # 7) Final label by priority
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
        # 8-phase logic based on cycle_fraction
        #  We define each phase in 1/8 increments
        if cycle_fraction < 0.0625 or cycle_fraction >= 0.9375:
            label = "new_moon"
        elif cycle_fraction < 0.1875:
            label = "waxing_crescent"
        elif cycle_fraction < 0.3125:
            label = "first_quarter"
        elif cycle_fraction < 0.4375:
            label = "waxing_gibbous"
        elif cycle_fraction < 0.5625:
            label = "full_moon"
        elif cycle_fraction < 0.6875:
            label = "waning_gibbous"
        elif cycle_fraction < 0.8125:
            label = "last_quarter"
        else:
            label = "waning_crescent"

    debug_lines.append(f"Final label: {label}")

    # Print debug info to stderr (helpful if using GitHub Actions)
    for dbg in debug_lines:
        print(f"[DEBUG] {dbg}", file=sys.stderr)

    return label

if __name__ == "__main__":
    print(get_moon_label(), end="")
