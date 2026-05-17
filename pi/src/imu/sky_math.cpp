#include "imu/sky_math.h"

#include <cmath>

namespace imu::sky {

namespace {
constexpr double TAU = 6.283185307179586476925286766559;
constexpr double DEG2RAD = 0.017453292519943295769;
}

double local_sidereal_time(double unix_time, double lon_east) {
    // GMST polynomial from Meeus, "Astronomical Algorithms", eq. 12.4.
    double JD  = unix_time / 86400.0 + 2440587.5;
    double D   = JD - 2451545.0;
    double T   = D / 36525.0;
    double gmst_deg = 280.46061837
                    + 360.98564736629 * D
                    + 0.000387933 * T * T
                    - T * T * T / 38710000.0;
    double lst_deg = gmst_deg + lon_east;
    double lst_rad = std::fmod(lst_deg * DEG2RAD, TAU);
    if (lst_rad < 0) lst_rad += TAU;
    return lst_rad;
}

void world_to_cel(float lat_north, float mag_decl, float lst, float R[9]) {
    // Build R_trueNED→cel first, then post-multiply by R_magNED→trueNED so
    // the resulting matrix maps the magnetometer-referenced world frame
    // directly into celestial.
    //
    // True-NED basis vectors expressed in celestial coords (φ = lat_north,
    // θ = lst):
    //   north_cel = [-sin φ cos θ, -sin φ sin θ,  cos φ]
    //   east_cel  = [-sin θ,        cos θ,        0  ]
    //   down_cel  = [-cos φ cos θ, -cos φ sin θ, -sin φ]
    float cp = std::cos(lat_north), sp = std::sin(lat_north);
    float ct = std::cos(lst),       st = std::sin(lst);

    float Tn[9] = {
        -sp * ct, -st, -cp * ct,
        -sp * st,  ct, -cp * st,
         cp,       0,  -sp
    };

    // R_magNED→trueNED: rotation around the down axis by +mag_decl (so that
    // a vector along magnetic north — body's x in the world frame — lands
    // `mag_decl` east of true north).
    float cd = std::cos(mag_decl), sd = std::sin(mag_decl);
    float D[9] = {
        cd, -sd, 0,
        sd,  cd, 0,
         0,   0, 1
    };

    matmul3(Tn, D, R);
}

void radec_roll_to_matrix(float ra, float dec, float roll, float R[9]) {
    float cr = std::cos(ra),  sr = std::sin(ra);
    float cd = std::cos(dec), sd = std::sin(dec);
    float cp = std::cos(roll), sp = std::sin(roll);

    // Pointing direction (camera +Z) in celestial frame.
    float zx = cd * cr, zy = cd * sr, zz = sd;
    // Celestial north tangent at (ra, dec).
    float nx = -sd * cr, ny = -sd * sr, nz = cd;
    // Celestial east tangent at (ra, dec) (direction of increasing RA).
    float ex = -sr, ey = cr, ez = 0.0f;

    // At roll=0: camera x = east, camera y = north. Positive roll rotates
    // the image counter-clockwise as seen looking out from the camera, which
    // means image-up swings toward image-left (i.e. away from east).
    //   y_cam_cel = -sin(φ) east + cos(φ) north
    //   x_cam_cel =  cos(φ) east + sin(φ) north
    float xx = cp * ex + sp * nx;
    float xy = cp * ey + sp * ny;
    float xz = cp * ez + sp * nz;
    float yx = -sp * ex + cp * nx;
    float yy = -sp * ey + cp * ny;
    float yz = -sp * ez + cp * nz;

    // Columns = camera basis vectors in celestial frame.
    R[0] = xx; R[1] = yx; R[2] = zx;
    R[3] = xy; R[4] = yy; R[5] = zy;
    R[6] = xz; R[7] = yz; R[8] = zz;
}

void matrix_to_radec_roll(const float R[9],
                          float &ra, float &dec, float &roll) {
    float zx = R[2], zy = R[5], zz = R[8];
    // Clamp to keep asin sane against tiny float drift.
    float zclamp = zz; if (zclamp > 1.0f) zclamp = 1.0f; else if (zclamp < -1.0f) zclamp = -1.0f;
    dec = std::asin(zclamp);
    ra  = std::atan2(zy, zx);
    if (ra < 0) ra += static_cast<float>(6.283185307179586);

    // Reconstruct east/north at the pointing direction and project the camera
    // image-up axis (column 1 of R) onto them to recover roll.
    float cr = std::cos(ra), sr = std::sin(ra);
    float cd = std::cos(dec), sd = std::sin(dec);
    float nx = -sd * cr, ny = -sd * sr, nz = cd;
    float ex = -sr,      ey =  cr,      ez = 0.0f;

    float ux = R[1], uy = R[4], uz = R[7];
    float north_c = ux * nx + uy * ny + uz * nz;  //  cos(roll)
    float east_c  = ux * ex + uy * ey + uz * ez;  // -sin(roll)
    roll = std::atan2(-east_c, north_c);
}

void matrix_to_radec(const float R[9], float &ra, float &dec) {
    float zz = R[8];
    if (zz > 1.0f) zz = 1.0f; else if (zz < -1.0f) zz = -1.0f;
    dec = std::asin(zz);
    ra  = std::atan2(R[5], R[2]);
    if (ra < 0) ra += static_cast<float>(6.283185307179586);
}

void quat_to_matrix(const float q[4], float R[9]) {
    float w = q[0], x = q[1], y = q[2], z = q[3];
    R[0] = 1 - 2*(y*y + z*z);  R[1] = 2*(x*y - w*z);      R[2] = 2*(x*z + w*y);
    R[3] = 2*(x*y + w*z);      R[4] = 1 - 2*(x*x + z*z);  R[5] = 2*(y*z - w*x);
    R[6] = 2*(x*z - w*y);      R[7] = 2*(y*z + w*x);      R[8] = 1 - 2*(x*x + y*y);
}

void matmul3(const float A[9], const float B[9], float C[9]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            C[i*3 + j] = A[i*3 + 0] * B[0*3 + j]
                       + A[i*3 + 1] * B[1*3 + j]
                       + A[i*3 + 2] * B[2*3 + j];
        }
    }
}

void transpose3(const float R[9], float R_out[9]) {
    R_out[0] = R[0]; R_out[1] = R[3]; R_out[2] = R[6];
    R_out[3] = R[1]; R_out[4] = R[4]; R_out[5] = R[7];
    R_out[6] = R[2]; R_out[7] = R[5]; R_out[8] = R[8];
}

}  // namespace imu::sky
