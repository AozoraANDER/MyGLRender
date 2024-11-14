#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <vector> 
#include <algorithm>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <limits>

#define M_PI 3.14159265358979323846

std::vector<Vec2f> uv_coords;

Vec3f cross(const Vec3f &v1, const Vec3f &v2) {
  return Vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
               v1.x * v2.y - v1.y * v2.x);
}

// 计算叉积函数
void loadUVCoords(const std::string &filename, std::vector<Vec2f> &uv_coords) {
  std::ifstream in(filename, std::ifstream::in);
  if (in.fail()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    char trash;
    if (!line.compare(0, 3, "vt ")) {
      iss >> trash >> trash;  // skip "vt"
      Vec2f uv;
      iss >> uv.x >> uv.y;
      uv_coords.push_back(uv);
    }
  }
}

std::vector<std::vector<int>> face_uv_indices;

void loadFaceUVIndices(const std::string &filename,
                       std::vector<std::vector<int>> &face_uv_indices) {
  std::ifstream in(filename, std::ifstream::in);
  if (in.fail()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream iss(line);
    char trash;
    if (!line.compare(0, 2, "f ")) {
      iss >> trash;  // skip 'f'
      std::vector<int> face_uv;
      int v_idx, t_idx, n_idx;
      for (int i = 0; i < 3; i++) {
        char separator;
        iss >> v_idx >> separator >> t_idx >> separator >> n_idx;
        face_uv.push_back(t_idx - 1);  // 索引从1开始，改为从0开始
      }
      face_uv_indices.push_back(face_uv);
    }
  }
}

// 计算重心坐标
Vec3f barycentric(Vec2i A, Vec2i B, Vec2i C, Vec2i P) {
    Vec3f s[2];
    s[0] = Vec3f(C.x - A.x, B.x - A.x, A.x - P.x);
    s[1] = Vec3f(C.y - A.y, B.y - A.y, A.y - P.y);
    Vec3f u = cross(s[0], s[1]);  // 使用叉积函数计算

    if (std::abs(u.z) < 1) return Vec3f(-1, 1, 1);  // 三角形退化的情况

    return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

// 绘制三角形，使用深度缓冲
void triangle(Vec3f *world_coords, Vec2i *screen_coords, Vec2f *uv_coords,
              float *zbuffer, TGAImage &image, TGAImage &diffuse_map,
              float intensity) {
  // 你的包围盒计算和深度缓冲代码保持不变

    Vec2f bboxmin(std::numeric_limits<float>::max(),
                std::numeric_limits<float>::max());
  Vec2f bboxmax(-std::numeric_limits<float>::max(),
                -std::numeric_limits<float>::max());
  Vec2f clamp(image.get_width() - 1, image.get_height() - 1);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      bboxmin[j] = std::max(
          0.f, std::min(bboxmin[j], static_cast<float>(screen_coords[i][j])));
      bboxmax[j] = std::min(
          clamp[j],
          std::max(bboxmax[j], static_cast<float>(screen_coords[i][j])));

    }
  }

  // 在三角形内采样纹理颜色
  

  // 计算 P 点的深度值 (Z 值)
  Vec3f P;
  for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
    for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
      // 使用barycentric重心坐标的部分
      Vec3f bc_screen = barycentric(screen_coords[0], screen_coords[1],
                                    screen_coords[2], Vec2i(P.x, P.y));
      if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;

      // 计算 P 点的深度值 (Z 值)
      P.z = 0;
      for (int i = 0; i < 3; i++) {
        P.z += world_coords[i].z * bc_screen[i];
      }

      // 更新zbuffer并设置像素颜色
      int idx = static_cast<int>(P.x) + static_cast<int>(P.y) * image.get_width();
      if (zbuffer[idx] < P.z) {
        zbuffer[idx] = P.z;

        // 计算纹理坐标
        Vec2f uv(0, 0);
        for (int i = 0; i < 3; i++) {
          uv.x += uv_coords[i].x * bc_screen[i];
          uv.y += uv_coords[i].y * bc_screen[i];
        }
        uv.x *= diffuse_map.get_width();
        uv.y *= diffuse_map.get_height();

        // 应用光照强度
        TGAColor color =
            diffuse_map.get(static_cast<int>(uv.x), static_cast<int>(uv.y));
        color.bgra[0] = (unsigned char)(color.bgra[0] * intensity);
        color.bgra[1] = (unsigned char)(color.bgra[1] * intensity);
        color.bgra[2] = (unsigned char)(color.bgra[2] * intensity);

       Vec3f bc_screen = barycentric(screen_coords[0], screen_coords[1],
                                      screen_coords[2], Vec2i(P.x, P.y));

        if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) {
         return;  // 退出当前函数
         }
        image.set(P.x, P.y, color);
      }
    }
  }
}

int main() {
  TGAImage image(800, 800, TGAImage::RGB);
  int width = 800;
  int height = 800;
  Model model("assets/african_head.obj");

  // 加载UV坐标和面UV索引
  loadUVCoords("assets/african_head.obj", uv_coords);
  loadFaceUVIndices("assets/african_head.obj", face_uv_indices);

  TGAImage diffuse_map;
  if (!diffuse_map.read_tga_file("assets/african_head_diffuse.tga")) {
    std::cerr << "Error: Cannot load the texture file!" << std::endl;
    return 1;
  }
  diffuse_map.flip_vertically();

  float *zbuffer = new float[width * height];
  for (int i = 0; i < width * height; ++i) {
    zbuffer[i] = std::numeric_limits<float>::lowest();  // 初始设为负无穷
  }

  Vec3f light_dir(0, 0, -1);  // 定义光照方向

  for (int i = 0; i < model.nfaces(); i++) {
    std::vector<int> face = model.face(i);
    std::vector<int> face_uv = face_uv_indices[i];  // 获取每个面的UV索引

    Vec3f world_coords[3];
    Vec2i screen_coords[3];
    Vec2f uv_coords_triangle[3];

    for (int j = 0; j < 3; j++) {
      Vec3f v = model.vert(face[j]);
      world_coords[j] = v;
      screen_coords[j] =
          Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
      uv_coords_triangle[j] = uv_coords[face_uv[j]];  // 获取每个顶点的UV坐标
    }

    Vec3f n = cross(world_coords[2] - world_coords[0],
                    world_coords[1] - world_coords[0]);

    n.normalize();
    float intensity = n * light_dir;  // 计算光照强度

    if (intensity > 0) {
      triangle(world_coords, screen_coords, uv_coords_triangle, zbuffer, image,
               diffuse_map, intensity);
     }
    }

  image.flip_vertically();  // 翻转图像，使原点在左下角
  image.write_tga_file("Alpha.tga");

  delete[] zbuffer;  // 释放 zbuffer
  return 0;
}
