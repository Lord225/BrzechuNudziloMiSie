#include <iostream>
#include <Windows.h>
#include <algorithm>
#include <chrono>
#include <math.h>
#include "glm/vec2.hpp"
#include "glm/vec3.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/matrix.hpp"
#include "glm/gtc/type_ptr.hpp"

#define uint unsigned int 
#define samplemethod 0
#define COLOR_PALLETTE " .:-=+*#%@" //" 0123456789"
#define ROTATION_SPEED 1.0f
#define LIGHT_ROTATION_SPEED 1.5f


class ASCIIDisplay {
public:
	
	ASCIIDisplay(int H = 100, int W = 100) : H(H), W(W)
	{
		membuffor = new char[(W+1) * H];
		Zbuffer = new float[W * H];
	}
	void imShow() 
	{
		using namespace std;
		cout << membuffor << endl;
	}
	void imDrawPoint(glm::vec3 pos, float val)
	{
		auto vert = view_to_screen_space(pos);
		setpixel(vert.x, vert.y, get_color(val), pos.z);
	}
	void imClear(float val) 
	{
		char val_char = get_color(val);
		for (uint i = 0; i < H; i++)
		{
			for (uint j = 0; j < W; j++)
			{
				membuffor[i*(W+1) + j] = val_char;
				Zbuffer[i * (W) + j] = -100000;
			}
			membuffor[i *(W+1) + W] = '\n';
		}
		membuffor[(W + 1) * H - 1] = '\0';
	}

	void imDrawTriangle(glm::vec3 vert[], glm::vec3 colorA, glm::vec3 colorB, glm::vec3 colorC, float (*shader)(glm::vec2, glm::vec3))
	{
		drawTriangle(vert, colorA, colorB, colorC, shader);
	}

	//arr (vertex array): (3x float (pos) + 3x float (color))xN , ind - indicates, (size)
	void imDraw(const float * arr, int start, int size, glm::mat4 transform, glm::mat3 inversedTransposedModel, float (*shader)(glm::vec2, glm::vec3))
	{
		int drawcalls = 0;
		glm::vec3 pos[3], color[3];
		for (int i = start; i < size; i+=3)
		{
			for (int j = 0; j < 3; j++)
			{
				pos[j][0] = arr[6 * (i + j) + 0];
				pos[j][1] = arr[6 * (i + j) + 1];
				pos[j][2] = arr[6 * (i + j) + 2];

				pos[j] = transform * glm::vec4(pos[j], 1.0f);
			

				color[j][0] = arr[6 * (i + j) + 3];
				color[j][1] = arr[6 * (i + j) + 4];
				color[j][2] = arr[6 * (i + j) + 5];

				color[j] = inversedTransposedModel * color[j];
			}	
			
			//pos[1] = { 0.0f, 1.0f, -0.2f };
			//pos[2] = { -0.5f, -0.75f, 1.5f };
			//pos[0] = { -0.5f, 0.75f, 0.5f };
			
			imDrawTriangle(pos, color[0], color[1], color[2], shader);

			//imDrawPixel(pos[0], 1.0);
			//imDrawPixel(pos[1], 1.0);
			//imDrawPixel(pos[2], 1.0);
			
			drawcalls++;
		}
		//std::cout << drawcalls << std::endl;
	}
private:
	void setpixel(int x, int y, float val, float z)
	{
		setpixel(x, y, get_color(val), z);
	}
	
	char get_color(float val)
	{
#if samplemethod == 0
		//Linear Clip
		return COLOR_MAP[(int)(std::clamp(val, 0.0f, 1.0f) * (strlen(COLOR_MAP) - 1))];

#elif samplemethod == 1
		//Sigmoid
		return COLOR_MAP[(int)(sigmoid(val) * (strlen(COLOR_MAP) - 1))];

#elif samplemethod == 2
		//Linear Soft Clip
		return COLOR_MAP[(int)(r_elu(val) * (strlen(COLOR_MAP) - 1))];
#endif
	}

	glm::vec3 view_to_screen_space(glm::vec3 pos)
	{
		pos = (glm::clamp(pos, { -1.0f, -1.0f, -1.0 }, { 1.0f, 1.0f, 1.0f }) + 1.0f) / 2.0f;
		pos.x = pos.x * (W - 1);
		pos.y = pos.y * (H - 1);
		return pos;
	}
	
	void drawTriangle(glm::vec3 vert[], glm::vec3 colorA, glm::vec3 colorB, glm::vec3 colorC, float (* shade)(glm::vec2, glm::vec3))
	{
		glm::vec2 bbmin = { INFINITY, INFINITY }, bbmax = { -INFINITY, -INFINITY };
		glm::vec3 vproj[3];
		for (int i = 0; i < 3; ++i) {
			vproj[i] = view_to_screen_space(vert[i]);
			vproj[i].z = vert[i].z;
			if (vproj[i].x < bbmin.x) bbmin.x = vproj[i].x;
			if (vproj[i].y < bbmin.y) bbmin.y = vproj[i].y;
			if (vproj[i].x > bbmax.x) bbmax.x = vproj[i].x;
			if (vproj[i].y > bbmax.y) bbmax.y = vproj[i].y;
		}

		uint xmin = std::max(0, std::min(W - 1, (int)std::floor(bbmin.x)));
		uint ymin = std::max(0, std::min(H - 1, (int)std::floor(bbmin.y)));
		uint xmax = std::max(0, std::min(W - 1, (int)std::floor(bbmax.x)));
		uint ymax = std::max(0, std::min(H - 1, (int)std::floor(bbmax.y)));

		float max = 1000;

		for (int y = ymin; y <= ymax; ++y) {
			for (int x = xmin; x <= xmax; ++x) {
				// check of if current pixel lies in triangle
				if (PointInTriangle(x, y, vproj[0], vproj[1], vproj[2])) {
					float v, w, u;
					get_barycentric(x, y, vproj[0], vproj[1], vproj[2], v, w, u);
					float depth = (vproj[0].z * v + vproj[1].z * w + vproj[2].z * u) / 3;
					auto color = (colorA * v + colorB * w + colorC * u)/3.0f;

					if (getZ(x, y) < depth) {
						setpixel(x, y, shade({ x, y }, color), depth);
					}
				}
			}
		}
		//std::cout << max << std::endl;
	}

	void get_barycentric(float x, float y, glm::vec3 A, glm::vec3 B, glm::vec3 C, float & v, float & w, float & u) {
		glm::vec2 v0 = B - A, v1 = C - A, v2 = glm::vec3(x, y, 0) - A;
		float d00 = dot(v0, v0);
		float d01 = dot(v0, v1);
		float d11 = dot(v1, v1);
		float d20 = dot(v2, v0);
		float d21 = dot(v2, v1);
		float denom = d00 * d11 - d01 * d01;
		v = (d11 * d20 - d01 * d21) / denom;
		w = (d00 * d21 - d01 * d20) / denom;
		u = 1.0f - v - w;
	}

	float area(float x1, float y1, float x2, float y2, float x3, float y3)
	{
		return abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0f);
	}

	bool PointInTriangle(float x, float y, glm::vec2 A, glm::vec2 B, glm::vec2 C)
	{
		float Aa = area(A.x, A.y, B.x, B.y, C.x, C.y);

		float A1 = area(x, y, B.x, B.y, C.x, C.y);

		float A2 = area(A.x, A.y, x, y, C.x, C.y);

		float A3 = area(A.x, A.y, B.x, B.y, x, y);

		return (abs(Aa-(A1+A2+A3)) < 0.001f);
	}

	void setpixel(int x, int y, char val, float z)
	{
		membuffor[y*(W+1)+ x] = val;
		Zbuffer[y * (W + 1) + x] = z;
	}
	float getZ(int x, int y) {
		return Zbuffer[y * (W + 1) + x];
	}

	float sigmoid(float x) {
		return 1 / (1 + exp(-5*x+3));
	}
	float r_elu(float x) {
		return x > 0 ? 1 - exp(-2 * x) : 0;
	}
	
	const char* COLOR_MAP = COLOR_PALLETTE;
	int H, W;
	char* membuffor;
	float* Zbuffer;
	
};

class PerspectiveCamera {
public:
	glm::vec3 pos;
	PerspectiveCamera(float fov = 110, float aspects = 1.0f) 
	{
		pos = { 0, 0, 3.0f };
		projection = glm::perspective(glm::radians(fov), aspects, 0.01f, 200.0f);
	}
	//Transforms model into view space
	glm::mat4 getTransformMatrix(glm::mat4 model)
	{
		view = glm::translate(glm::mat4(1.0f), pos);
		return projection * view * model;
	}
	//For normal calculations
	glm::mat3 getInversedTransposedModel(glm::mat4 model)
	{
		return glm::mat3(glm::transpose(glm::inverse(model)));
	}
private:
	glm::mat4 projection;
	glm::mat4 view;
};

//deltaTime and total time tracker
class TimeMessurmet {
public:
	
	TimeMessurmet(float fps_limit = 16.6f){
		last_time = std::chrono::system_clock::now();
		delta = 1.0;
		total_time = 0;
		this->fps_limit = fps_limit;
	}

	void update() {
		

		if (fps_limit != 0) {
			Sleep(std::max(0.0, fps_limit - delta));
		}
		now = std::chrono::system_clock::now();
		delta = std::chrono::duration_cast<std::chrono::microseconds>((now - last_time)).count()/1'000'000.0f;
		last_time = now;
		total_time += delta;
	}

	float deltaTime() {
		return delta;
	}
	double deltaTime_d() {
		return delta;
	}
	double time() {
		return total_time;
	}

private:
	double delta;
	double total_time;
	float fps_limit;
	std::chrono::system_clock::time_point last_time;
	std::chrono::system_clock::time_point now;
};

glm::vec3 light_pos = {0.0f, 0.0, 0.0};


//Shader used to draw cube
float shader(glm::vec2 pos, glm::vec3 color) 
{
	return glm::clamp(glm::dot(color, light_pos)+0.4f, 0.2f, 1.0f);
}

glm::mat4 model(1.0f);

int main()
{
	auto engine = ASCIIDisplay(30, 80);
	auto camera = PerspectiveCamera();
	auto time = TimeMessurmet();
	
	// x,y, z, normal.x, normal.y, normal.z
	const float vertices[] = {
	-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
	 0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
	0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,

	0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
	-0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,
	-0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,

	 0.5f, -0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
	-0.5f, -0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
	 0.5f,  0.5f,  0.5f,  0.0f,  0.0f, 1.0f,

	 0.5f,  0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
	-0.5f, -0.5f,  0.5f,  0.0f,  0.0f, 1.0f,
	-0.5f,  0.5f,  0.5f,  0.0f,  0.0f, 1.0f,

	-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
	-0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
	-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
	-0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,
	-0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,
	-0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,

	 0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
	 0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
	 0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
	 0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,
	 0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,
	 0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,

	-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
	 0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,
	 0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
	 0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
	-0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,
	-0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,

	-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
	 0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,
	 0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
	 0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
	-0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,
	-0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f
	};
	
	while (true) {
		engine.imClear(0.0f);
		
		auto transform = camera.getTransformMatrix(model);
		auto inversed = camera.getInversedTransposedModel(model);

		model = glm::rotate(model, ROTATION_SPEED * time.deltaTime(), { 1.0f,0.5f,0.0f });

		//draw vertex array
		engine.imDraw(vertices, 0, 36, transform, inversed, shader);
		engine.imDrawPoint(light_pos, 1.0f);

		system("cls");

		//move light in circle
		light_pos.x = sin(time.time() * LIGHT_ROTATION_SPEED) * 1;
		light_pos.y = cos(time.time() * LIGHT_ROTATION_SPEED) * 1;
		

		engine.imShow();
		time.update();
		//std::cout << time.deltaTime() << std::endl;

	}
}
