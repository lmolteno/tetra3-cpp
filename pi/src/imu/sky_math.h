#pragma once

// Coordinate-system glue for the IMU-driven pointing estimate.
//
// Three frames matter:
//   - Body / camera frame:  z = optical axis out, y = image up, x = image right
//   - World (magnetic NED): x = magnetic north, y = east, z = down
//   - Celestial (equatorial): x = vernal equinox, z = north celestial pole
//
// We treat the IMU body frame and camera frame as identical for the initial
// pre-solve estimate; after the first plate solve a fixed R_cam→imu correction
// is learned and applied.
//
// Rotation matrices are stored row-major 3x3 floats (R[0..8]).

namespace imu::sky {

// Local sidereal time, radians in [0, 2π).
// unix_time : seconds since 1970-01-01T00:00:00Z (decimal allowed)
// lon_east  : observer longitude in degrees, positive east
double local_sidereal_time(double unix_time, double lon_east);

// Build R_world→celestial(t).
// lat_north  : observer latitude  (radians, positive north)
// mag_decl   : magnetic declination (radians, positive east — i.e. magnetic
//              north lies `mag_decl` east of true north). Composed with the
//              true-NED↔celestial rotation so the world frame is the actual
//              IMU/magnetometer frame (mag-aligned).
// lst        : local sidereal time (radians, see local_sidereal_time)
void world_to_cel(float lat_north, float mag_decl, float lst, float R[9]);

// Build R_cam→celestial from (RA, Dec, roll). Convention: camera frame has
// z = optical axis, y = image up, x = image right. At roll=0 the image's up
// direction aligns with celestial north and right aligns with celestial east
// (i.e. increasing RA). Positive roll rotates the image counter-clockwise as
// seen from behind the camera.
void radec_roll_to_matrix(float ra, float dec, float roll, float R[9]);

// Inverse of radec_roll_to_matrix.
void matrix_to_radec_roll(const float R[9],
                          float &ra, float &dec, float &roll);

// Extract just the pointing direction (RA, Dec) — equivalent to
// matrix_to_radec_roll but skipping the roll work.
void matrix_to_radec(const float R[9], float &ra, float &dec);

// Quaternion [w,x,y,z] → 3x3 rotation matrix s.t. v_world = R * v_body
// (Madgwick / Hamilton convention).
void quat_to_matrix(const float q[4], float R[9]);

// C = A * B (row-major 3x3 matrix multiply).
void matmul3(const float A[9], const float B[9], float C[9]);

// R_out = R^T (transpose, equivalently inverse for an orthogonal R).
void transpose3(const float R[9], float R_out[9]);

}  // namespace imu::sky
