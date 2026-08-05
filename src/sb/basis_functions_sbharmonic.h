/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    ***       Karamelo       ***
 *               Parallel Material Point Method Simulator
 * 
 * Copyright (2019) Alban de Vaucorbeil, alban.devaucorbeil@monash.edu
 * Materials Science and Engineering, Monash University
 * Clayton VIC 3800, Australia

 * This software is distributed under the GNU General Public License.
 *
 * ----------------------------------------------------------------------- */

#ifndef MPM_BASIS_FUNCTIONS_H
#define MPM_BASIS_FUNCTIONS_H

#include "D:\code\code-gitee\myLib\funclib.h"
using namespace std;

namespace BasisFunction { 
    // ---------- 2维坐标转换:笛卡尔->局部比例边界 ----------
    inline bool rectanglePhysicalToSBLocal(
        const double x, const double y, const double xc, const double yc,
        const double dx, const double dy,
        int& localEdgeIdx, double& xi, double& eta) noexcept //承诺不会抛出任何异常
    {
        constexpr double tol = 1e-12;
        // 归一化坐标u,v -> [-1,1]
        const double u = 2.0 * (x - xc) / dx;
        const double v = 2.0 * (y - yc) / dy;
        const double au = std::abs(u);
        const double av = std::abs(v);
        // 计算径向坐标,判断是否超出单元,防止÷0
        xi = au >= av ? au : av;
        if (xi > 1.0 + tol) return false;
        if (xi < tol) {
            localEdgeIdx = 0;
            xi = 0.0;
            eta = 0.0;
            return true;
        }
        // 计算局部边索引和周向坐标
        if (au >= av) {// ? R : L
            localEdgeIdx = u >= 0.0 ? 1 : 3;
            eta = v / u;
        }
        else {// ? B : T
            localEdgeIdx = v >= 0.0 ? 2 : 0;
            eta = -u / v;
        }
        return true;
    }

	// ---------- 3维坐标转换:笛卡尔->局部比例边界 ----------
	enum class AbaqusHexFace : int {
		S1 = 1, // z-，节点 1-2-3-4，内法向 +z
		S2 = 2, // z+，节点 5-8-7-6，内法向 -z
		S3 = 3, // y-，节点 1-5-6-2，内法向 +y
		S4 = 4, // x+，节点 2-6-7-3，内法向 -x
		S5 = 5, // y+，节点 3-7-8-4，内法向 -y
		S6 = 6  // x-，节点 4-8-5-1，内法向 +x
	};
	/**
	 * 笛卡尔坐标转换为长方体比例边界坐标。
	 * 假设：
	 *   1. 单元为轴对齐长方体；
	 *   2. 缩放中心为 (xc, yc, zc)；
	 *   3. dx、dy、dz 为长方体完整边长；
	 *   4. 面编号采用 Abaqus C3D8 的 S1～S6；
	 *   5. eta-zeta 参数面的法向指向单元内部。
	 * 输出：
	 *   localFaceIdx = Abaqus 面编号 1～6
	 *   xi           = 径向坐标，0 为缩放中心，1 为边界
	 *   eta, zeta    = 面内坐标，范围 [-1,1]
	 */
	inline bool cuboidPhysicalToSBLocal(
		const double x, const double y, const double z,
		const double xc, const double yc, const double zc,
		const double dx, const double dy, const double dz,
		int& localFaceIdx, double& xi, double& eta, double& zeta) noexcept
	{
		constexpr double tol = 1.0e-12;
		// 归一化笛卡尔坐标, u,v,w ∈ [-1,1]
		const double u = 2.0 * (x - xc) / dx;
		const double v = 2.0 * (y - yc) / dy;
		const double w = 2.0 * (z - zc) / dz;
		const double au = std::abs(u);
		const double av = std::abs(v);
		const double aw = std::abs(w);
		// L∞ 径向坐标,判断是否超出单元,防止÷0
		xi = std::max({ au, av, aw });
		if (xi > 1.0 + tol) { return false; }
		// 缩放中心处不存在唯一的所属面,指定为第一个面 S1.
		if (xi < tol) {
			localFaceIdx = static_cast<int>(AbaqusHexFace::S1);
			xi = 0.0;
			eta = 0.0;
			zeta = 0.0;
			return true;
		}
		/*
		* 面内参数方向经过选择，使：
		*     d(x_b)/d(eta) × d(x_b)/d(zeta)
		* 指向单元内部。
		* 边和角点同时属于多个面，因此必须规定并列优先级。
		* 当前优先级为：
		*     x 面 > y 面 > z 面
		* 即与原二维实现中 au >= av 的处理方式一致。
		*/
		if (au >= av && au >= aw) { // x 主方向：S4 或 S6
			if (u >= 0.0) {// Abaqus S4：x = dx/2, 内法向：-x
				localFaceIdx = static_cast<int>(AbaqusHexFace::S4);
				eta = -v / u;
				zeta = w / u;
			}
			else {// Abaqus S6：x = -dx/2, 内法向：+x
				localFaceIdx = static_cast<int>(AbaqusHexFace::S6);
				eta = -v / u;
				zeta = -w / u;
			}
		}
		else if (av >= aw) {// y 主方向：S3 或 S5
			if (v >= 0.0) {// Abaqus S5：y = +dy/2, 内法向：-y
				localFaceIdx = static_cast<int>(AbaqusHexFace::S5);
				eta = u / v;
				zeta = w / v;
			}
			else {// Abaqus S3：y = -dy/2,内法向：+y
				localFaceIdx = static_cast<int>(AbaqusHexFace::S3);
				eta = -u / v;
				zeta = -w / v;
			}
		}
		else {// z 主方向：S1 或 S2
			if (w >= 0.0) {// Abaqus S2：z = +dz/2,内法向：-z
				localFaceIdx = static_cast<int>(AbaqusHexFace::S2);
				eta = u / w;
				zeta = -v / w;
			}
			else {// Abaqus S1：z = -dz/2, 内法向：+z
				localFaceIdx = static_cast<int>(AbaqusHexFace::S1);
				eta = -u / w;
				zeta = -v / w;
			}
		}
		// 处理浮点舍入造成的轻微越界。
		xi = std::clamp(xi, 0.0, 1.0);
		eta = std::clamp(eta, -1.0, 1.0);
		zeta = std::clamp(zeta, -1.0, 1.0);
		return true;
	}

	// ---------- 计算形函数值 ----------

	// ----------  ----------

	// ----------  ----------

	// ----------  ----------

	// ----------  ----------

}
#endif
