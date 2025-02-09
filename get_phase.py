import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import get_sun, get_moon, GeocentricTrueEcliptic
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
      7) else one of the basic 8 phases:
         'new_moon', 'waxing_crescent', 'first_quarter', 'waxing_gibbous',
         'full_moon', 'waning_gibbous', 'last_quarter', 'waning_crescent'.

    This is a demo. Real astronomy for eclipses, equinox-based harvest moons,
    etc. can be more complex.
    """

    # Current UTC time
    now = Time.now()

    # Geocentric positions of Sun and Moon
    moon_pos = get_moon(now)
    sun_pos = get_sun(now)

    # 1) PHASE FRACTION: (1 + cos angular separation) / 2
    elongation = moon_pos.separation(sun_pos)
    phase_fraction = (1 + np.cos(elongation.to(u.rad))) / 2  # 0..1

    # 2) DISTANCE to the Moon (km)
    distance_km = moon_pos.distance.to(u.km).value

    # 3) ECLIPTIC LAT (for rough eclipse check)
    moon_ecliptic = moon_pos.transform_to(GeocentricTrueEcliptic(equinox=now))
    ecliptic_lat_deg = moon_ecliptic.lat.to(u.deg).value

    # Define thresholds
    FULL_MOON_THRESHOLD = 0.95   # fraction >= 0.95 => near full
    NEW_MOON_THRESHOLD  = 0.05   # fraction <= 0.05 => near new
    SUPERMOON_PERIGEE   = 360000 # ~360,000 km as a rough cutoff for "supermoon"
    MICRO_MOON_APOGEE   = 405000 # ~405,000 km as a rough cutoff for "micro_moon"
    ECLIPTIC_LAT_THRESH = 1.5    # deg from ecliptic => possible eclipse if near full

    is_near_full  = (phase_fraction > FULL_MOON_THRESHOLD)
    is_near_new   = (phase_fraction < NEW_MOON_THRESHOLD)
    is_near_eclip = (abs(ecliptic_lat_deg) < ECLIPTIC_LAT_THRESH)

    # 4) Rough Eclipse Check (blood moon)
    #    For a real eclipse check, you'd compute Earth’s umbra geometry.
    #    We'll just say "lunar_eclipse" if near full + near ecliptic.
    might_be_eclipse = is_near_full and is_near_eclip

    # 5) Supermoon / Micro Moon (if near full)
    #    Common definitions vary. We'll say "supermoon" if < 360,000 km
    #    from Earth AND near full; "micro_moon" if > 405,000 km AND near full.
    is_supermoon  = (distance_km < SUPERMOON_PERIGEE and is_near_full)
    is_micro_moon = (distance_km > MICRO_MOON_APOGEE  and is_near_full)

    # 6) Harvest Moon (near full in late Sept or early Oct)
    d_utc = now.to_datetime()  # Python datetime in UTC
    is_harvest_time = (d_utc.month == 9 or d_utc.month == 10)
    is_harvest_moon = (is_harvest_time and is_near_full)

    # 7) Blue Moon (second full moon in the same calendar month)
    #    We'll do a quick backward search to see if there's another near-full
    #    within this month. Approximate, but works OK in practice.
    is_blue_moon = False
    if is_near_full:
        current_month = d_utc.month
        day_step = 0.5  # half-day increments
        t_search = now
        found_another_full = False
        for _ in range(60):  # 60 half-day steps => up to 30 days
            t_search = t_search - day_step*u.day
            dt_utc = t_search.to_datetime()
            if dt_utc.month != current_month:
                # We stepped into a previous month
                break
            # Check if near-full at this step
            mpos = get_moon(t_search)
            spos = get_sun(t_search)
            el = mpos.separation(spos)
            frac = (1 + np.cos(el.to(u.rad))) / 2
            if frac > FULL_MOON_THRESHOLD:
                found_another_full = True
                break
        is_blue_moon = found_another_full

    # 8) Black Moon (second new moon in the same calendar month)
    #    Similar approach: if is_near_new, look backward to see if there's
    #    another near-new in the same month.
    is_black_moon = False
    if is_near_new:
        current_month = d_utc.month
        day_step = 0.5
        t_search = now
        found_another_new = False
        for _ in range(60):  # up to ~30 days
            t_search = t_search - day_step*u.day
            dt_utc = t_search.to_datetime()
            if dt_utc.month != current_month:
                break
            # Check if near-new at this step
            mpos = get_moon(t_search)
            spos = get_sun(t_search)
            el = mpos.separation(spos)
            frac = (1 + np.cos(el.to(u.rad))) / 2
            if frac < NEW_MOON_THRESHOLD:
                found_another_new = True
                break
        is_black_moon = found_another_new

    # 9) Priority or Combining Labels
    #    A single date can be "Super Blue Blood Moon," etc.
    #    Here we pick ONLY ONE label in a chain:
    if might_be_eclipse:
        return "lunar_eclipse"       # a.k.a. "blood_moon"
    elif is_blue_moon:
        return "blue_moon"
    elif is_black_moon:
        return "black_moon"
    elif is_harvest_moon:
        return "harvest_moon"
    elif is_supermoon:
        return "supermoon"
    elif is_micro_moon:
        return "micro_moon"

    # Otherwise, pick from the 8 standard phases:
    if is_near_new:
        return "new_moon"
    elif phase_fraction < 0.1875:
        return "waxing_crescent"
    elif phase_fraction < 0.3125:
        return "first_quarter"
    elif phase_fraction < 0.4375:
        return "waxing_gibbous"
    elif phase_fraction < 0.5625:
        return "full_moon"
    elif phase_fraction < 0.6875:
        return "waning_gibbous"
    elif phase_fraction < 0.8125:
        return "last_quarter"
    elif phase_fraction < 0.9375:
        return "waning_crescent"
    else:
        return "new_moon"  # near 1.0 again

if __name__ == "__main__":
    print(get_moon_label(), end="")
