#ifndef CHAOSMAP_GENERATOR_HPP_
#define CHAOSMAP_GENERATOR_HPP_

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "vector3d.hpp"

class ChaosMapGenerator {
   public:
    ChaosMapGenerator() {}
    virtual ~ChaosMapGenerator() {}

    virtual void frame_transformation() = 0;
    virtual void trajectory_propagation(double timestep) = 0;
    virtual void calc_SALImap_params(double delta_V_lowerlimit,
                                     double delta_V_upperlimit,
                                     double delta_V_changerate,
                                     double timestep,
                                     double time_end) = 0;
};

class CRTBPChaosMap : public ChaosMapGenerator {
   private:
    //! init_ast_j2000_[0]:慣性座標系における小惑星の初期位置init_ast_kanseie_state_[1]:慣性座標系における小惑星の初期位置での速度
    std::array<Vector3d, 2> init_e_G_;
    //! init_ast_j2000_[0]:慣性座標系における小惑星の初期位置init_ast_kanseie_state_[1]:慣性座標系における小惑星の初期位置での速度
    std::array<Vector3d, 2> init_ast_G_;
    //! init_ast_rotframe_state_:回転座標系における小惑星の初期位置
    std::array<Vector3d, 2> init_ast_R_state_;
    //! trajectory_params_:軌道のパラメータ; [0]:t, [1]:x, [2]:y, [3]:z, [4]:vx, [5]:vy,
    //! [6]:vz
    std::vector<std::array<double, 7>> trajectory_params_;
    //! chaos_map_params_:カオスマップのパラメータ; [0]:t, [1]:x, [2]:y, [3]:z, [4]:delta V,
    //! [5]~[19]:SALI
    std::vector<std::array<double, 20>> SALI_map_params_;
    double SoI_radius_;
    double FA_radius_;
    double mu_;

    static constexpr double G = 6.67430e-11;
    static constexpr double AU = 1.49597870700e11;
    static constexpr double M_SUN = 1.989400e30;
   public:
    /**
     * @brief こんとらくた
     * @param init_e_j2000 j2000系における地球の初期位置
     * @param init_ast_j2000_ j2000系における小惑星の初期位置
     * @param SoI_radius the radius of Sphere of Influence
     * @param FA_radius the radius of forbidden area
     * @param mu 重力定数
     */


    CRTBPChaosMap(std::array<Vector3d, 2> init_e_j2000,
                  std::array<Vector3d, 2> init_ast_j2000,
                  double SoI_radius = 0.03,
                  double FA_radius = 0.00007,
                  double mu = 3.003e-6);
    /**
     * @brief
     * 慣性座標系での小惑星の初期位置init_ast_j2000_を
     * 太陽-地球重心の回転座標系から見た初期位置init_ast_rotframe_state_へ変換
     * @return void
     */
    void frame_transformation() override;
    void frame_transformation_();

    /**
     * @brief crtbpを用いてSoI_radiusを抜けるまで軌道計算
     * @return void
     */
    void trajectory_propagation(double timestep) override;

    // /**
    //  * @brief 任意の点において速さを変更する
    //  * @param delta_V 速さの変化量
    //  */
    // void modify_ast_speed(double delta_V) override;

    /**
     * @brief trajectory_propagationで計算した軌道上の各点においてSALIを計算
     * chaos_map_params_の末尾にSALIを追加
     * @return void
     */
    void calc_SALImap_params(double delta_V_lowerlimit,
                             double delta_V_upperlimit,
                             double delta_V_changerate,
                             double timestep,
                             double time_end) override;

    void write_trajectory_params2file(const std::string& filename);
    void write_SALImap_params2file(const std::string& filename);
};

#endif  // CHAOSMAP_GENERATOR_HPP_