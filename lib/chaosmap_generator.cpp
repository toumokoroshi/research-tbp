#include "chaosmap_generator.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "crtbp.hpp"
#include "vector3d.hpp"

#define _debug
static std::array<std::array<double, 3>, 3> create_frame_conversion_matrix(const Vector3d& unit_n,
                                                                           double theta) {
    double half_theta = theta / 2.0;
    double q0 = std::cos(half_theta);
    double sin_half_theta = std::sin(half_theta);
    double q1 = unit_n.x() * sin_half_theta;
    double q2 = unit_n.y() * sin_half_theta;
    double q3 = unit_n.z() * sin_half_theta;

    double q0q0 = q0 * q0;
    double q1q1 = q1 * q1;
    double q2q2 = q2 * q2;
    double q3q3 = q3 * q3;
    double q0q1 = q0 * q1;
    double q0q2 = q0 * q2;
    double q0q3 = q0 * q3;
    double q1q2 = q1 * q2;
    double q1q3 = q1 * q3;
    double q2q3 = q2 * q3;

    // 回転行列
    std::array<std::array<double, 3>, 3> convert_matrix = {
        {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 + q0q3), 2.0 * (q1q3 - q0q2)},
         {2.0 * (q1q2 - q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 + q0q1)},
         {2.0 * (q1q3 + q0q2), 2.0 * (q2q3 - q0q1), q0q0 - q1q1 - q2q2 + q3q3}}
    };
    return convert_matrix;
};

// 回転軸と回転角からクオータニオン経由で回転行列を生成
static std::array<std::array<double, 3>, 3> create_coordinate_rotation_matrix(
    const Vector3d& unit_n, double theta) {
    double half_theta = theta / 2.0;
    double q0 = std::cos(half_theta);
    double sin_half_theta = std::sin(half_theta);
    double q1 = unit_n.x() * sin_half_theta;
    double q2 = unit_n.y() * sin_half_theta;
    double q3 = unit_n.z() * sin_half_theta;

    double q0q0 = q0 * q0;
    double q1q1 = q1 * q1;
    double q2q2 = q2 * q2;
    double q3q3 = q3 * q3;
    double q0q1 = q0 * q1;
    double q0q2 = q0 * q2;
    double q0q3 = q0 * q3;
    double q1q2 = q1 * q2;
    double q1q3 = q1 * q3;
    double q2q3 = q2 * q3;

    std::array<std::array<double, 3>, 3> rot_matrix = {
        {{q0q0 + q1q1 - q2q2 - q3q3, 2.0 * (q1q2 - q0q3), 2.0 * (q1q3 + q0q2)},
         {2.0 * (q1q2 + q0q3), q0q0 - q1q1 + q2q2 - q3q3, 2.0 * (q2q3 - q0q1)},
         {2.0 * (q1q3 - q0q2), 2.0 * (q2q3 + q0q1), q0q0 - q1q1 - q2q2 + q3q3}}
    };
    return rot_matrix;
};

// ロドリゲス回転公式を表す行列を生成するラムダ関数
static std::array<std::array<double, 3>, 3> create_rodrigues_matrix(const Vector3d& unit_n,
                                                                    double theta) {
    double c = std::cos(theta);
    double s = std::sin(theta);
    double t = 1.0 - c;
    double x = unit_n.x();
    double y = unit_n.y();
    double z = unit_n.z();

    std::array<std::array<double, 3>, 3> rot_matrix = {
        {{t * x * x + c, t * x * y - s * z, t * x * z + s * y},
         {t * x * y + s * z, t * y * y + c, t * y * z - s * x},
         {t * x * z - s * y, t * y * z + s * x, t * z * z + c}}
    };
    return rot_matrix;
};

// ベクトルに回転行列を適用するラムダ関数
static Vector3d apply_matrix(std::array<std::array<double, 3>, 3> rot_matrix, const Vector3d& v) {
    return {rot_matrix[0][0] * v.x() + rot_matrix[0][1] * v.y() + rot_matrix[0][2] * v.z(),
            rot_matrix[1][0] * v.x() + rot_matrix[1][1] * v.y() + rot_matrix[1][2] * v.z(),
            rot_matrix[2][0] * v.x() + rot_matrix[2][1] * v.y() + rot_matrix[2][2] * v.z()};
};

CRTBPChaosMap::CRTBPChaosMap(std::array<Vector3d, 2> init_e_G,
                             std::array<Vector3d, 2> init_ast_G,
                             double SoI_radius,
                             double FA_radius,
                             double mu)
    : init_e_G_(init_e_G),
      init_ast_G_(init_ast_G),
      SoI_radius_(SoI_radius),
      FA_radius_(FA_radius),
      mu_(mu) {}

void CRTBPChaosMap::frame_transformation() {
    /* j2000系から見ているとき :_G */

    // クオータニオンの微分wotukatte kaitenngyouretuno bibun wo keisan
    auto Cdot = [](std::pair<double, Vector3d> q,
                   std::pair<double, Vector3d> qdot) -> std::array<std::array<double, 3>, 3> {
        double q0 = q.first;
        Vector3d qv = q.second;
        double qdot0 = qdot.first;
        Vector3d qdotv = qdot.second;

        std::array<std::array<double, 3>, 3> Cdot_ = {
            {{2.0 * (q0 * qdot0 + qv.x() * qdotv.x() - qv.y() * qdotv.y() - qv.z() * qdotv.z()),
              2.0 * (qv.x() * qdotv.y() + qdotv.x() * qv.y() + q0 * qdotv.z() + qdot0 * qv.z()),
              2.0 * (qv.x() * qdotv.z() + qdotv.x() * qv.z() - q0 * qdotv.y() - qdot0 * qv.y())},
             {2.0 * (qv.x() * qdotv.y() + qdotv.x() * qv.y() - q0 * qdotv.z() - qdot0 * qv.z()),
              2.0 * (q0 * qdot0 - qv.x() * qdotv.x() + qv.y() * qdotv.y() - qv.z() * qdotv.z()),
              2.0 * (qv.y() * qdotv.z() + qdotv.y() * qv.z() + q0 * qdotv.x() + qdot0 * qv.x())},
             {2.0 * (qv.x() * qdotv.z() + qdotv.x() * qv.z() + q0 * qdotv.y() + qdot0 * qv.y()),
              2.0 * (qv.y() * qdotv.z() + qdotv.y() * qv.z() - q0 * qdotv.x() - qdot0 * qv.x()),
              2.0 * (q0 * qdot0 - qv.x() * qdotv.x() - qv.y() * qdotv.y() + qv.z() * qdotv.z())}}
        };
        return Cdot_;
    };

    // 入力されたベクトルを入力された行列で変換
    auto convert = [](std::array<std::array<double, 3>, 3> convert_matrix,
                      const Vector3d& v) -> Vector3d {
        return {convert_matrix[0][0] * v.x() + convert_matrix[0][1] * v.y() +
                    convert_matrix[0][2] * v.z(),
                convert_matrix[1][0] * v.x() + convert_matrix[1][1] * v.y() +
                    convert_matrix[1][2] * v.z(),
                convert_matrix[2][0] * v.x() + convert_matrix[2][1] * v.y() +
                    convert_matrix[2][2] * v.z()};
    };

    Vector3d G_x{1., 0., 0.};
    Vector3d G_y{0., 1., 0.};
    Vector3d G_z{0., 0., 1.};
    /* 慣性系のz軸を地球の軌道面の法線ベクトルに一致させるための回転軸を含む面を計算*/
    // 地球の軌道面の法線ベクトル
    Vector3d norm_n_e_G = init_e_G_[0].gaiseki(init_e_G_[1]).normalise();
    //  慣性系のz軸を地球の軌道面の法線ベクトルに一致させるための回転軸を含む面の法線ベクトル
    Vector3d n1{norm_n_e_G.x(), norm_n_e_G.y(), norm_n_e_G.z() - 1.};

    /* 慣性系のx軸を慣性系での地球ベクトルに一致させるための回転軸を含む面を計算*/
    Vector3d norm_init_e_G_pos = init_e_G_[0].normalise();
    Vector3d n2{norm_init_e_G_pos.x() - 1., norm_init_e_G_pos.y(), norm_init_e_G_pos.z()};
    /* 回転軸を決定 */
    Vector3d rotax_buf_G = n1.gaiseki(n2).normalise();
    /*  回転角の決定と回転マトリックスの計算*/
    std::array<std::array<double, 3>, 3> convert_G_to_R = {0.0};
    std::array<std::array<double, 3>, 3> convert_R_to_kansei = {0.0};
    std::array<std::array<double, 3>, 3> Cdot_G_to_R = {0.0};

    // Vector3d p = norm_n_e_G - rotax_buf * (rotax_buf.naiseki(norm_n_e_G));
    // rotaxが示す面に投影されたnorm_n_e_G
    //  Vector3d q = G_z - rotax_buf * (rotax_buf.naiseki(G_z));
    // rotaxを法線ベクトルに持つ面に投影されたnorm_init_e_G_pos
    Vector3d p = norm_init_e_G_pos - rotax_buf_G * (rotax_buf_G.naiseki(norm_init_e_G_pos));
    // rotaxが示す面に投影されたG_x
    Vector3d q = G_x - rotax_buf_G * (rotax_buf_G.naiseki(G_x));
    //    rotaxが示す面に投影されたG_z
    double theta = std::atan2(p.gaiseki(q).magnitude(), p.naiseki(q));

    std::array<double, 3> buf = {rotax_buf_G.x(), rotax_buf_G.y(), rotax_buf_G.z()};
    if (rotax_buf_G.gaiseki(q).naiseki(p) < 0) {
        theta = -theta;
    }
    Vector3d rotax_G{buf[0], buf[1], buf[2]};
    convert_G_to_R = create_frame_conversion_matrix(rotax_G, theta);
    convert_R_to_kansei = create_coordinate_rotation_matrix(rotax_G, theta);

    /* 速度変換の準備 */
    double G = 6.67430e-11;
    double AU = 1.49597870700e11;
    double M_SUN = 1.989400e30;
    double ND_time_ref = std::sqrt(AU * AU * AU / (G * M_SUN));
    double mean_motion = 1 / ND_time_ref;
    // AU/day -> 無次元
    // Vector3d ND_init_e_velo_G = init_e_G_[1] * ND_time_ref / 86400.0;
    Vector3d ND_init_e_velo_G = norm_n_e_G.gaiseki(init_e_G_[0]).normalise() *
                                init_e_G_[1].magnitude() * ND_time_ref / 86400.;
    double po = init_e_G_[0].magnitude();
    Vector3d dot_norm_r_e =
        (ND_init_e_velo_G * po * po - init_e_G_[0] * (ND_init_e_velo_G.naiseki(init_e_G_[0]))) /
        po / po / po;
    double npx = norm_n_e_G.x();
    double npy = norm_n_e_G.y();
    double npz = norm_n_e_G.z();
    double ex = norm_init_e_G_pos.x();
    double ey = norm_init_e_G_pos.y();
    double ez = norm_init_e_G_pos.z();
    double vex = dot_norm_r_e.x();
    double vey = dot_norm_r_e.y();
    double vez = dot_norm_r_e.z();

    double n1crossn2_norm = n1.gaiseki(n2).magnitude();
    Vector3d dot_n1crossn2{
        npy * vez - (npz - 1.) * vey, (npz - 1.) * vex - npx * vez, npx * vey - npy * vex};
    double dot_n1crossn2norm = ((npy * npy * ez * vez + (1. - npz) * (1. - npz) * ey * vey +
                                 npy * (1. - npz) * (vez * ey + ez * vey)) +
                                (1. - npz) * (1. - npz) * (ex - .1) * vex + npx * npx * ez * vez -
                                (1. - npz) * npx * (-vex * ez + (1. - ex) * vez) +
                                (npx * npx * ey * vey - npy * npy * (1. - ex) * vex +
                                 npx * npy * (vey * (1. - ex) - ey * vex))) /
                               n1crossn2_norm;

    Vector3d dot_rotax = (dot_n1crossn2 * n1crossn2_norm - n1.gaiseki(n2) * dot_n1crossn2norm) /
                         n1crossn2_norm / n1crossn2_norm;

    double dot_pdotq = vex + (rotax_G.naiseki(norm_init_e_G_pos) * dot_rotax.x()) +
                       (dot_norm_r_e.naiseki(rotax_G) + norm_init_e_G_pos.naiseki(dot_rotax));
    double dot_pnorm =
        (norm_init_e_G_pos.naiseki(dot_norm_r_e) +
         norm_init_e_G_pos.naiseki(rotax_G) *
             (dot_norm_r_e.naiseki(rotax_G) + norm_init_e_G_pos.naiseki(dot_rotax))) /
        p.magnitude();
    double dot_qnorm = rotax_G.x() * dot_rotax.x() / q.magnitude();

    double dot_theta = -(dot_pdotq * p.magnitude() * q.magnitude() -
                         p.naiseki(q) * (dot_pnorm * q.magnitude() + p.magnitude() * dot_qnorm)) /
                       p.magnitude() / p.magnitude() / q.magnitude() / q.magnitude() /
                       std::sin(theta);
    Vector3d q_buf =
        rotax_G * dot_theta * std::cos(theta / 2.) / 2. + dot_rotax * std::sin(theta / 2.);
    std::pair<double, Vector3d> qdot = {-dot_theta * std::sin(theta / 2.) / 2., q_buf};
    std::pair<double, Vector3d> q__ = {std::cos(theta / 2.), rotax_G * std::sin(theta / 2.)};

    Cdot_G_to_R = Cdot(q__, qdot);

    /* 小惑星の位置の座標変換 */

    Vector3d e_R_pos = convert(convert_G_to_R, init_e_G_[0]);
    Vector3d e_R_velo =
        convert(Cdot_G_to_R, init_e_G_[0]) + convert(convert_G_to_R, ND_init_e_velo_G);
    // Vector3d e_R_velo = convert(Cdot_G_to_R, init_e_G_[0]) +
    //                       convert(convert_G_to_R,init_e_G_[1]);
    Vector3d ast_R_pos_buf = convert(convert_G_to_R, init_ast_G_[0]);
    Vector3d ast_velo_buf = convert(Cdot_G_to_R, init_ast_G_[0]) +
                            convert(convert_G_to_R, init_ast_G_[1] * ND_time_ref / 86400.0);

    // Vector3d G_z_R = convert(convert_R_to_kansei, G_z);
    Vector3d G_z_R = convert(convert_R_to_kansei, G_z);
    // //  check
    std::cout << std::setprecision(15) << "converted earth position = " << e_R_pos.x() << " "
              << e_R_pos.y() << " " << e_R_pos.z() << std::endl;
    std::cout << std::setprecision(15) << "converted z ax of j2000 = " << G_z_R.x() << " "
              << G_z_R.y() << " " << G_z_R.z() << std::endl;
    std::cout << std::setprecision(15) << "normal vector of earth orbit j2000= " << norm_n_e_G.x()
              << " " << norm_n_e_G.y() << " " << norm_n_e_G.z() << std::endl;
    // ast_R_pos_buf
    std::cout << std::setprecision(15) << "converted asteroid position = " << ast_R_pos_buf.x()
              << " " << ast_R_pos_buf.y() << " " << ast_R_pos_buf.z() << std::endl;

    std::cout << std::setprecision(15) << "initial earth velocity = " << ND_init_e_velo_G.x() << " "
              << ND_init_e_velo_G.y() << " " << ND_init_e_velo_G.z() << std::endl;
    std::cout << std::setprecision(15) << "converted earth velocity = " << e_R_velo.x() << " "
              << e_R_velo.y() << " " << e_R_velo.z() << std::endl;
    std::cout << std::setprecision(15) << "converted asteroid velocity = " << ast_velo_buf.x()
              << " " << ast_velo_buf.y() << " " << ast_velo_buf.z() << std::endl;

    // 地球の位置がx＝1になるように調整
    Vector3d ast_R_pos{
        ast_R_pos_buf.x() + 1.0 - e_R_pos.x() - mu_, ast_R_pos_buf.y(), ast_R_pos_buf.z()};

    init_ast_R_state_[0] = ast_R_pos;
    init_ast_R_state_[1] = ast_velo_buf;
    std::cout << "<>               Asteroid data converted to rotating frame" << std::endl;
}

void CRTBPChaosMap::frame_transformation_() {
    Vector3d G_x{1.0, 0.0, 0.0};
    Vector3d G_y{0.0, 1.0, 0.0};
    Vector3d G_z{0.0, 0.0, 1.0};
    /*x軸を一致させる中間座標系rot1と移す  */
    // 慣性系から見た地球ベクトルと慣性系のx軸がなす角度を計算
    double theta1 = std::atan2(init_e_G_[0].gaiseki(G_x).magnitude(), init_e_G_[0].naiseki(G_x));
    std::cout << std::setprecision(15);

    // 回転軸を決定
    Vector3d rotax1 = G_x.gaiseki(init_e_G_[0]).normalise();
    // 変換行列を生成
    std::array<std::array<double, 3>, 3> convert_G_to_R1 = create_rodrigues_matrix(rotax1, -theta1);
    std::array<std::array<double, 3>, 3> rotate1 = create_rodrigues_matrix(rotax1, theta1);

    /* 地球の軌道面の法線ベクトルと中間座標系rot1のz軸を一致させ目的の回転座標系rot2へ移す */
    // 地球の軌道面の法線ベクトル
    Vector3d n_plane = init_e_G_[0].gaiseki(init_e_G_[1]).normalise();

    // 慣性系から見た中間座標系rot1のz軸
    Vector3d rot1_z = apply_matrix(rotate1, G_z);

    //    慣性座標系から見た地球軌道面の法線ベクトルと慣性系から見た中間座標系rot1のz軸がなす角度を計算
    double theta2 = std::atan2(n_plane.gaiseki(rot1_z).magnitude(), n_plane.naiseki(rot1_z));
    Vector3d rotax2 = {1.0, 0.0, 0.0};

    // 変換行列を生成
    std::array<std::array<double, 3>, 3> convert_R1_to_R2 =
        create_rodrigues_matrix(rotax2, -theta2);
    std::array<std::array<double, 3>, 3> rotate2 = create_rodrigues_matrix(rotax2, theta2);
    // 小惑星の位置を回転座標系に変換
    Vector3d ast_pos_R1 = apply_matrix(convert_G_to_R1, init_ast_G_[0]);
    Vector3d epos_R1 = apply_matrix(convert_G_to_R1, init_e_G_[0]);
    Vector3d ast_pos_R2 = apply_matrix(convert_R1_to_R2, ast_pos_R1);
    Vector3d epos_R2 = apply_matrix(convert_R1_to_R2, epos_R1);
    Vector3d ast_pos_R{ast_pos_R2.x() + 1.0 - epos_R2.x() - mu_, ast_pos_R2.y(), ast_pos_R2.z()};

    /* calc the velocity of rotating frame */
    //  変換行列convert_G_to_R1の微分を計算
    // theta1の微分

    double ND_time_ref = std::sqrt(AU * AU * AU / (G * M_SUN));  // Non-Dimensional time
    double mean_motion = 1 / ND_time_ref;
    std::cout << "mean motion = " << mean_motion << std::endl;
    // conbert AU/day to Non-Dimensional
    Vector3d ND_init_e_velo_G =
        n_plane.gaiseki(init_e_G_[0]).normalise() * init_e_G_[1].magnitude() * ND_time_ref / 86400.;

    double e_G_x = init_e_G_[0].x();
    double e_G_y = init_e_G_[0].y();
    double e_G_z = init_e_G_[0].z();
    double ND_e_G_vx = ND_init_e_velo_G.x();
    double ND_e_G_vy = ND_init_e_velo_G.y();
    double ND_e_G_vz = ND_init_e_velo_G.z();

    double poo = init_e_G_[0].magnitude();
    double sin_theta1 = std::sqrt(1. - std::cos(theta1) * std::cos(theta1));

    // double dot_theta1 = -(ND_e_G_vx * poo * poo -
    //                       e_G_x * (e_G_x * ND_e_G_vx + e_G_y * ND_e_G_vy + e_G_z * ND_e_G_vz)) /
    //                       <- 円制限三体問題だからちきゅうの速度と位置の内積は0
    //                     poo / poo / poo / sin_theta1;
    double dot_theta1 = -ND_e_G_vx / poo / sin_theta1;

    // theta2の微分
    double sin_theta2 = std::sqrt(1. - std::cos(theta2) * std::cos(theta2));
    if (n_plane.gaiseki(apply_matrix(rotate1, G_x)).naiseki(rot1_z) < 0) {
        sin_theta2 = -sin_theta2;
    }

    Vector3d k{n_plane.x(), n_plane.y(), n_plane.z()};

    double n1 = std::sqrt(e_G_y*e_G_y +e_G_z*e_G_z);
    double a = e_G_y * ND_e_G_vy + e_G_z * ND_e_G_vz;
    double b = ND_e_G_vy * e_G_z + e_G_y * ND_e_G_vz;

    double dot_theta2 = (k.x() * (ND_e_G_vz * n1 * n1 - e_G_z * a) * sin_theta1 / n1 / n1 / n1 -
                         k.x() * e_G_z * std ::cos(theta1) * dot_theta1 / n1 +
                         k.y() * (b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1. - std::cos(theta1)) /
                             n1 / n1 / n1 / n1 +
                         k.y() * e_G_z * ND_e_G_vz * sin_theta1 * dot_theta1 / n1 / n1 -
                         2. * k.z() * (e_G_y * ND_e_G_vy * n1 * n1 - e_G_y * a) *
                             (1. - std::cos(theta1)) / n1 / n1 / n1 / n1 -
                         k.z() * e_G_y * e_G_y * sin_theta1 * dot_theta1 / n1 / n1 +
                         k.z() * sin_theta1 * dot_theta1) /
                        sin_theta2;

    double c1 = G_x.naiseki(init_e_G_[0]) / init_e_G_[0].magnitude();
    double s1 = sin_theta1;
    double c2 = std::cos(theta2);
    double s2 = sin_theta2;
    std::array<std::array<double, 3>, 3> dot_convert_G_to_R1{
        {{-s1 * dot_theta1,
          (ND_e_G_vy * n1 * n1 - e_G_y * a) * s1 / n1 / n1 / n1 + e_G_y * c1 * dot_theta1 / n1,
          (ND_e_G_vz * n1 * n1 - e_G_z * a) * s1 / n1 / n1 / n1 + e_G_z * c1 * dot_theta1 / n1},
         {-(ND_e_G_vy * n1 * n1 - e_G_y * a) * s1 / n1 / n1 / n1 - e_G_y * c1 * dot_theta1 / n1,
          (2. * e_G_z * ND_e_G_vz * n1 * n1 - 2. * e_G_z * e_G_z * a) * (1.0 - c1) / n1 / n1 / n1 /
                  n1 +
              e_G_z * e_G_z * s1 * dot_theta1 / n1 / n1 - s1 * dot_theta1,
          -(b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1. - c1) / n1 / n1 / n1 / n1 -
              e_G_z * e_G_y * s1 * dot_theta1 / n1 / n1},
         {-(ND_e_G_vz * n1 * n1 - e_G_z * a) * s1 / n1 / n1 / n1 - e_G_z * c1 * dot_theta1 / n1,
          -(b * n1 * n1 - 2. * e_G_z * e_G_y * a) * (1.0 - c1) / n1 / n1 / n1 / n1 -
              e_G_z * e_G_y * s1 * dot_theta1 / n1 / n1,
          (2. * e_G_y * ND_e_G_vy * n1 * n1 - e_G_y * e_G_y * a) * (1. - c1) / n1 / n1 / n1 / n1 +
              e_G_y * e_G_y * s1 * dot_theta1 / n1 / n1 - s1 * dot_theta1}}
    };

    std::array<std::array<double, 3>, 3> dot_convert_R1_to_R2{
        {{0.0, 0.0, 0.0},
         {0.0, -s2 * dot_theta2, c2 * dot_theta2},
         {0.0, -c2 * dot_theta2, -s2 * dot_theta2}}
    };
    // 速度を座標変換
    Vector3d ast_vel_R_1 =
        apply_matrix(dot_convert_R1_to_R2, apply_matrix(convert_G_to_R1, init_ast_G_[0]));
    Vector3d ast_vel_R_2 =
        apply_matrix(convert_R1_to_R2, apply_matrix(dot_convert_G_to_R1, init_ast_G_[0]));
    Vector3d ast_vel_R_3 =
        apply_matrix(convert_R1_to_R2, apply_matrix(convert_G_to_R1, init_ast_G_[1]));

    Vector3d ast_vel_R = ast_vel_R_1 + ast_vel_R_2 + ast_vel_R_3;
    // check
    Vector3d e_vel_R_1 = {0, 0, 0};
    // apply_matrix(dot_convert_R1_to_R2, apply_matrix(convert_G_to_R1, init_e_G_[0]));
    Vector3d e_vel_R_2 = {0, 0, 0};
    // apply_matrix(convert_R1_to_R2, apply_matrix(dot_convert_G_to_R1, init_e_G_[0]));
    Vector3d e_vel_R_3 =
        apply_matrix(convert_R1_to_R2, apply_matrix(convert_G_to_R1, ND_init_e_velo_G));
    Vector3d e_vel_R = e_vel_R_1 + e_vel_R_2 + e_vel_R_3;

#ifdef _debug

    // 中間座標系にから見た地球の速度
    // OK!!!!!

#endif

    // check
    std::cout << "converted earth velocity = " << e_vel_R.x() << " " << e_vel_R.y() << " "
              << e_vel_R.z() << std::endl;
    init_ast_R_state_ = {ast_pos_R, ast_vel_R};

    std::cout << "<>               Asteroid data converted to rotating frame" << std::endl;
}

void CRTBPChaosMap::trajectory_propagation(double timestep) {
    // CRTBPクラスのインスタンスを生成
    CRTBP crtbp(init_ast_R_state_, timestep);
    double time = 0.0;
    std::cout << "<>sdfsd = " << crtbp.calc_r1() << std::endl;
    // 影響球を抜けるまで軌道計算
    while (crtbp.calc_r2() < SoI_radius_) {
        crtbp.symplectic_integration_step();
        time += timestep;
        // 軌道のパラメータを保存
        std::array<double, 7> params = {time,
                                        crtbp.current_state()[0],
                                        crtbp.current_state()[1],
                                        crtbp.current_state()[2],
                                        crtbp.current_state()[3],
                                        crtbp.current_state()[4],
                                        crtbp.current_state()[5]};
        trajectory_params_.push_back(params);
    }
}

void CRTBPChaosMap::calc_SALImap_params(double delta_V_lowerlimit,
                                        double delta_V_upperlimit,
                                        double delta_V_changerate,
                                        double timestep,
                                        double time_end) {
    // 任意の軌道要素の速度をdeltaV だけ変更するラムダ関数
    auto modify_ast_speed = [](const std::array<double, 7>& state,
                               double delta_V) -> std::array<double, 6> {
        std::array<double, 6> modified_state{
            state[1], state[2], state[3], state[4], state[5], state[6]};
        double speed_before =
            std::sqrt(state[4] * state[4] + state[5] * state[5] + state[6] * state[6]);
        modified_state[4] = modified_state[4] * (speed_before + delta_V) / speed_before;
        modified_state[5] = modified_state[5] * (speed_before + delta_V) / speed_before;
        modified_state[6] = modified_state[6] * (speed_before + delta_V) / speed_before;
        return modified_state;
    };

    // 任意の2組の軌道要素からSALIを計算するラムダ関数
    auto calc_SALI_from2pair = [](const std::array<double, 6>& ref_state,
                                  const std::array<double, 6>& perturbed_state1,
                                  const std::array<double, 6>& perturbed_state2) -> double {
        double SALI = 0.0;
        Vector3d deviation_vec1{perturbed_state1[0] - ref_state[0],
                                perturbed_state1[1] - ref_state[1],
                                perturbed_state1[2] - ref_state[2]};
        Vector3d deviation_vec2{perturbed_state2[0] - ref_state[0],
                                perturbed_state2[1] - ref_state[1],
                                perturbed_state2[2] - ref_state[2]};
        Vector3d normalise_dev_vec1 = deviation_vec1.normalise();
        Vector3d normalise_dev_vec2 = deviation_vec2.normalise();

        Vector3d sa = normalise_dev_vec1 - normalise_dev_vec2;
        Vector3d wa = normalise_dev_vec1 + normalise_dev_vec2;

        double sa_norm = sa.magnitude();
        double wa_norm = wa.magnitude();

        SALI = (sa_norm > wa_norm) ? wa_norm : sa_norm;
        return SALI;
    };

    /* 与える変更量を変えながらすべての軌道上の点に対してSALIを計算 */
    for (double delta_V = delta_V_lowerlimit; delta_V < delta_V_upperlimit;
         delta_V += delta_V_changerate) {
        /* 各軌道要素に対してSALIを計算 */
        for (auto& params : trajectory_params_) {
            // 速度を変更
            std::array<double, 6> modified_state = modify_ast_speed(params, delta_V);
            // 各パラメータに微小な擾乱を加える
            double delta = 1.0e-10;
            std::array<double, 6> x_perturbed_state{modified_state[0] + delta,
                                                    modified_state[1],
                                                    modified_state[2],
                                                    modified_state[3],
                                                    modified_state[4],
                                                    modified_state[5]};
            std::array<double, 6> y_perturbed_state{modified_state[0],
                                                    modified_state[1] + delta,
                                                    modified_state[2],
                                                    modified_state[3],
                                                    modified_state[4],
                                                    modified_state[5]};
            std::array<double, 6> z_perturbed_state{modified_state[0],
                                                    modified_state[1],
                                                    modified_state[2] + delta,
                                                    modified_state[3],
                                                    modified_state[4],
                                                    modified_state[5]};
            std::array<double, 6> vx_perturbed_state{modified_state[0],
                                                     modified_state[1],
                                                     modified_state[2],
                                                     modified_state[3] + delta,
                                                     modified_state[4],
                                                     modified_state[5]};
            std::array<double, 6> vy_perturbed_state{modified_state[0],
                                                     modified_state[1],
                                                     modified_state[2],
                                                     modified_state[3],
                                                     modified_state[4] + delta,
                                                     modified_state[5]};
            std::array<double, 6> vz_perturbed_state{modified_state[0],
                                                     modified_state[1],
                                                     modified_state[2],
                                                     modified_state[3],
                                                     modified_state[4],
                                                     modified_state[5] + delta};

            CRTBP crtbp(modified_state, timestep);
            CRTBP crtbp_x(x_perturbed_state, timestep);
            CRTBP crtbp_y(y_perturbed_state, timestep);
            CRTBP crtbp_z(z_perturbed_state, timestep);
            CRTBP crtbp_vx(vx_perturbed_state, timestep);
            CRTBP crtbp_vy(vy_perturbed_state, timestep);
            CRTBP crtbp_vz(vz_perturbed_state, timestep);

            bool abort_SALI = false;
            bool calc_SALI = true;
            double time = 0.0;
            // それぞれの軌道要素に対してCRTBPclassを生成し微小な擾乱を加えた軌道を計算
            while (calc_SALI) {
                if (time > time_end) {
                    break;
                } else {
                    crtbp.symplectic_integration_step();
                    crtbp_x.symplectic_integration_step();
                    crtbp_y.symplectic_integration_step();
                    crtbp_z.symplectic_integration_step();
                    crtbp_vx.symplectic_integration_step();
                    crtbp_vy.symplectic_integration_step();
                    crtbp_vz.symplectic_integration_step();

                    if (crtbp.calc_r1() < FA_radius_ || crtbp.calc_r1() > SoI_radius_) {
                        abort_SALI = true;
                        calc_SALI = false;
                        break;
                    } else {
                        time += timestep;
                    }
                }
            }

            std::array<double, 20> chaos_map_param;
            chaos_map_param[0] = params[0];
            chaos_map_param[1] = params[1];
            chaos_map_param[2] = params[2];
            chaos_map_param[3] = params[3];
            chaos_map_param[4] = delta_V;
            for (size_t i = 5; i < 20; i++) {
                chaos_map_param[i] = -1.0;
            }

            // 地球から離れすぎるか近すぎた場合は計算を中断、SALIを計算しない（デフォルト値－1のまま）
            if (abort_SALI) {
                break;
            } else {
                // 6つの微小な擾乱を加えた軌道要素から2つを選ぶ組み合わせ計15通りのSALIを計算
                std::array<std::array<double, 6>, 7> latest_perturbed_states{
                    crtbp.current_state(),
                    crtbp_x.current_state(),
                    crtbp_y.current_state(),
                    crtbp_z.current_state(),
                    crtbp_vx.current_state(),
                    crtbp_vy.current_state(),
                    crtbp_vz.current_state()};

                for (int i = 0; i < 6; i++) {
                    double index = 5;  // カオスマップのパラメータ; [0]:t, [1]:x, [2]:y,
                                       // [3]:z, [4]:delta Vなのでそれ以降を埋める
                    for (int j = i + 1; j < 6; j++) {
                        double SALI = calc_SALI_from2pair(latest_perturbed_states[0],
                                                          latest_perturbed_states[i],
                                                          latest_perturbed_states[j]);

                        chaos_map_param[index] = SALI;
                        index++;
                    }
                }
                SALI_map_params_.push_back(chaos_map_param);
            }
        }
    }
}

void CRTBPChaosMap::write_trajectory_params2file(const std::string& filename) {
    std::ofstream ofs(filename);
    for (auto& params : trajectory_params_) {
        for (auto& param : params) {
            ofs << param << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}

void CRTBPChaosMap::write_SALImap_params2file(const std::string& filename) {
    std::ofstream ofs(filename);
    for (auto& params : SALI_map_params_) {
        for (auto& param : params) {
            ofs << param << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}
