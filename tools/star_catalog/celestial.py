"""Quaternion-based orientation to celestial coordinate conversion.

Converts IMU quaternions (phone DeviceOrientation or BNO055) to RA/Dec/Roll
using astropy for sidereal time calculations.
"""

import math
from datetime import datetime, timezone

import numpy as np


def quaternion_to_rotation_matrix(q):
    """Convert quaternion [w, x, y, z] to 3x3 rotation matrix.

    World frame: X=East, Y=North, Z=Up.
    """
    w, x, y, z = q
    # Normalize
    n = math.sqrt(w * w + x * x + y * y + z * z)
    if n < 1e-10:
        return np.eye(3)
    w, x, y, z = w / n, x / n, y / n, z / n

    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


def quaternion_to_alt_az_roll(q):
    """Convert quaternion [w,x,y,z] to (altitude_deg, azimuth_deg, roll_deg).

    World frame: X=East, Y=North, Z=Up.
    Device frame: -Z = pointing direction (back of phone / camera).
    """
    R = quaternion_to_rotation_matrix(q)

    # Device -Z axis in world frame = pointing direction
    pointing = -R[:, 2]  # third column, negated
    px, py, pz = pointing

    # Altitude: angle above horizon
    alt_deg = math.degrees(math.asin(np.clip(pz, -1.0, 1.0)))

    # Azimuth: angle from North (Y) clockwise through East (X)
    az_deg = math.degrees(math.atan2(px, py)) % 360.0

    # Roll: rotation of device "up" around the pointing axis.
    # Device +Y axis in world frame = device "up" direction.
    device_up = R[:, 1]
    # Project device_up onto plane perpendicular to pointing direction
    # "expected up" = Z cross pointing, normalized, gives "right" in the sky plane
    # Then up_sky = pointing cross right
    horiz = np.array([-py, px, 0.0])  # pointing cross Z (simplified since Z=[0,0,1])
    horiz_len = np.linalg.norm(horiz)
    if horiz_len < 1e-6:
        # Pointing straight up or down — roll is ambiguous
        roll_deg = 0.0
    else:
        horiz /= horiz_len
        up_sky = np.cross(pointing, horiz)
        # Roll = angle between device_up projected and up_sky
        roll_deg = math.degrees(math.atan2(
            np.dot(device_up, horiz),
            np.dot(device_up, up_sky),
        ))

    return alt_deg, az_deg, roll_deg


def alt_az_to_ra_dec(alt_deg, az_deg, lat_deg, lon_deg, utc_time=None):
    """Convert horizontal coordinates to equatorial (RA, Dec) in degrees.

    Uses astropy for Local Sidereal Time calculation.
    """
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, AltAz, SkyCoord
    import astropy.units as u

    if utc_time is None:
        utc_time = datetime.now(timezone.utc)

    location = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg, height=0 * u.m)
    obstime = Time(utc_time)
    altaz_frame = AltAz(obstime=obstime, location=location)

    coord = SkyCoord(alt=alt_deg * u.deg, az=az_deg * u.deg, frame=altaz_frame)
    icrs = coord.icrs

    return icrs.ra.deg, icrs.dec.deg


def celestial_roll(q, lat_deg):
    """Compute roll relative to celestial north pole.

    Returns the angle (degrees) between the device's "up" direction and the
    direction toward the NCP, both projected onto the image plane
    (perpendicular to pointing direction).
    """
    R = quaternion_to_rotation_matrix(q)
    pointing = -R[:, 2]
    device_up = R[:, 1]

    # NCP unit vector in world frame (X=East, Y=North, Z=Up)
    lat_rad = math.radians(lat_deg)
    ncp = np.array([0.0, math.cos(lat_rad), math.sin(lat_rad)])

    # Project NCP onto image plane (perpendicular to pointing)
    ncp_proj = ncp - np.dot(ncp, pointing) * pointing
    ncp_len = np.linalg.norm(ncp_proj)
    if ncp_len < 1e-6:
        return 0.0  # pointing at/away from NCP
    ncp_proj /= ncp_len

    # Project device_up onto image plane
    up_proj = device_up - np.dot(device_up, pointing) * pointing
    up_len = np.linalg.norm(up_proj)
    if up_len < 1e-6:
        return 0.0
    up_proj /= up_len

    # Signed angle: NCP direction to device up, around pointing axis.
    # Negated to match the star map renderer convention where positive roll
    # rotates the star field so that NCP moves clockwise on screen.
    cross = np.cross(ncp_proj, up_proj)
    sin_angle = np.dot(cross, pointing)
    cos_angle = np.dot(ncp_proj, up_proj)

    return -math.degrees(math.atan2(sin_angle, cos_angle))


def imu_to_celestial(q, lat_deg, lon_deg, utc_time=None):
    """Full pipeline: quaternion + GPS + time -> (ra_rad, dec_rad, roll_rad).

    Roll is relative to celestial north (not local zenith).
    Returns radians to match pi_tracker's expected input format.
    """
    alt_deg, az_deg, _ = quaternion_to_alt_az_roll(q)
    ra_deg, dec_deg = alt_az_to_ra_dec(alt_deg, az_deg, lat_deg, lon_deg, utc_time)
    roll_deg = celestial_roll(q, lat_deg)

    return (
        math.radians(ra_deg),
        math.radians(dec_deg),
        math.radians(roll_deg),
    )
