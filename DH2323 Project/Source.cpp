#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "vector"
#include <algorithm>
#include <glm/glm.hpp>

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::vec2;
using glm::mat3;


// GLOBAL VARIABLES
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
const double f = SCREEN_HEIGHT * 0.99;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
float zBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

int t;
const int dataResX = 120; // size of height grid
const int dataResY = 120;
const int polyResX = (dataResX - 1) * 2; // number of polygons generated from height grid
const int polyResY = dataResY - 1;
const float scale = 0.25; // scale dictates size of polygons, how "big" the landscape will be

double h_persistence = 0.75; // noise generation parameters
int h_octave = 2;
double f_persistence = 0.75;
int f_octave = 2;
double t_persistence = 0.75;
int t_octave = 1;

vec3 cameraPos(0, 0, -3.001);
mat3 R(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
float yaw = 0;
float turnRate = 0.5;
vec3 lightPower = 1500.f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.1f * vec3(1, 1, 1);
vec3 lightPos((dataResX / 2)* scale, -5, (dataResY / 2)* scale);
int drawdistance = 25; // distance culling
int dimmingRange = 10; // dimming interval
SDL_Surface* screen;
float Lerp(float t, float a1, float a2);
vec3 Lerp(float t, vec3 a1, vec3 a2);
struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 pos3d, color, normal;
};

class Polygon {
public:
	vec3 v0, v1, v2, color, normal, center, lightv0, lightv1, lightv2, colorv0, colorv1, colorv2;

	float fertility, temperature;
	Polygon() {

	}
	Polygon(vec3 v0, vec3 v1, vec3 v2, float f, float t)
		:v0(v0), v1(v1), v2(v2)
	{
		center = (v0 + v1 + v2);
		center = vec3(center[0] / 3, center[1] / 3, center[2] / 3);
		fertility = f;
		temperature = t;
		ComputeNormal();

		setLights();

		vec3 colorUpper, colorLower;
		if (fertility < 0.20)
		{
			if (temperature > 0.6)
			{
				colorUpper = vec3(1, 0.6470, 0);//Color.Orange;
				colorLower = vec3(1, 0.2705, 0);//Color.OrangeRed;
			}
			else if (temperature > 0.4)
			{
				colorUpper = vec3(1, 0.6470, 0);//Color.Orange;
				colorLower = Lerp(0.25, vec3(0.5450, 0.2705, 0.0745), vec3(1, 0.2705, 0));
			}
			else
			{
				colorUpper = vec3(0.4, 0.8039, 0.6666);//Color.MediumAquamarine;
				colorLower = Lerp(0.5, vec3(0, 0.5450, 0.5450), vec3(0.4392, 0.5019, 0.5647)); //Color.Lerp(Color.DarkCyan, Color.SlateGray, 0.5f);
			}
		}
		else if (fertility < 0.40)
		{
			if (temperature > 0.8)
			{
				colorUpper = vec3(1, 0.6470, 0);//Color.Orange;
				colorLower = vec3(1, 0.2705, 0);//Color.OrangeRed;
			}
			else if (temperature > 0.6)
			{
				colorUpper = Lerp(0.45, vec3(0.5019, 0.6470, 0.1254), vec3(0.4862, 0.9882, 0));//Color.Lerp(Color.Goldenrod, Color.LawnGreen, 0.45f);
				colorLower = Lerp(0.35, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.4862, 0.9882, 0));
				//Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.LawnGreen, 0.35f);
			}
			else if (temperature > 0.4)
			{
				colorUpper = Lerp(0.75, vec3(0.5019, 0.6470, 0.1254), vec3(0.4196, 0.5568, 0.1372));//Color.Lerp(Color.Goldenrod, Color.OliveDrab, 0.75f);
				colorLower = Lerp(0.35, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.4196, 0.5568, 0.1372));//Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.OliveDrab, 0.35f);
			}
			else
			{
				colorUpper = Lerp(0.35, vec3(0.5019, 0.5019, 0), vec3(0.5960, 0.9843, 0.5960));//Color.Lerp(Color.Olive, Color.PaleGreen, 0.35f);
				colorLower = Lerp(0.25, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.5960, 0.9843, 0.5960));//Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.PaleGreen, 0.25f);
			}
		}
		else if (fertility < 0.60)
		{
			if (temperature > 0.6)
			{
				colorUpper = Lerp(0.35,
					Lerp(0.45, vec3(0.5019, 0.6470, 0.1254), vec3(0.4862, 0.9882, 0)),
					vec3(0.1333, 0.5450, 0.1333));//Color.Lerp(Color.Lerp(Color.Goldenrod, Color.LawnGreen, 0.45f), Color.ForestGreen, 0.35f);
				colorLower = Lerp(0.25,
					Lerp(0.15, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.4862, 0.9882, 0)),
					vec3(0.1333, 0.5450, 0.1333));//Color.Lerp(Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.LawnGreen, 0.15f), Color.ForestGreen, 0.25f);
			}
			else if (temperature > 0.4)
			{
				colorUpper = Lerp(0.3,
					Lerp(0.75, vec3(0.5019, 0.6470, 0.1254), vec3(0.4196, 0.5568, 0.1372)),
					vec3(0, 0.5019, 0));// Color.Lerp(Color.Lerp(Color.Goldenrod, Color.OliveDrab, 0.75f), Color.Green, 0.30f);
				colorLower = Lerp(0.15,
					Lerp(0.35, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.4196, 0.5568, 0.1372)),
					vec3(0, 0.5019, 0));// Color.Lerp(Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.OliveDrab, 0.35f), Color.Green, 0.15f);
			}
			else
			{
				colorUpper = Lerp(0.35,
					Lerp(0.35, vec3(0.5019, 0.5019, 0), vec3(0.5960, 0.9843, 0.5960)),
					vec3(0.1803, 0.5450, 0.3411));//Color.Lerp(Color.Lerp(Color.Olive, Color.PaleGreen, 0.35f), Color.SeaGreen, 0.35f);
				colorLower = Lerp(0.25,
					Lerp(0.15, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.5960, 0.9843, 0.5960)),
					vec3(0.1803, 0.5450, 0.3411));//Color.Lerp(Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.PaleGreen, 0.15f), Color.SeaGreen, 0.25f);
			}

		}
		else
		{
			if (temperature > 0.6)
			{
				colorUpper = Lerp(0.5,
					Lerp(0.5, vec3(0.5019, 0.6470, 0.1254), vec3(0.4862, 0.9882, 0)),
					vec3(0, 0.5019, 0));//Color.Lerp(Color.Lerp(Color.Goldenrod, Color.LawnGreen, 0.5f), Color.Green, 0.5f);
				colorLower = Lerp(0.35,
					Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)),
					vec3(0, 0.5019, 0));// Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.Green, 0.35f);
			}
			else if (temperature > 0.4)
			{
				colorUpper = Lerp(0.65,
					Lerp(0.35, vec3(0.5019, 0.6470, 0.1254), vec3(0.4196, 0.5568, 0.1372)),
					vec3(0, 0.5019, 0));//Color.Lerp(Color.Lerp(Color.Goldenrod, Color.OliveDrab, 0.35f), Color.Green, 0.65f);
				colorLower = Lerp(0.25,
					Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)),
					vec3(0.4196, 0.5568, 0.1372));// Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.OliveDrab, 0.25f);
			}
			else
			{
				colorUpper = Lerp(0.65,
					Lerp(0.35, vec3(0.5019, 0.5019, 0), vec3(0.5960, 0.9843, 0.5960)),
					vec3(0.1803, 0.5450, 0.3411));//Color.Lerp(Color.Lerp(Color.Olive, Color.PaleGreen, 0.35f), Color.SeaGreen, 0.65f);
				colorLower = Lerp(0.25,
					Lerp(0.15, Lerp(0.35, vec3(0, 0.3921, 0), vec3(0, 0, 0)), vec3(0.5960, 0.9843, 0.5960)),
					vec3(0, 0.5450, 0.5450));//Color.Lerp(Color.Lerp(Color.DarkGreen, Color.Black, 0.35f), Color.DarkCyan, 0.25f);
			}

		}
		colorv0 = Lerp(-v0.y * FastInverse(25), colorLower, colorUpper);
		colorv1 = Lerp(-v1.y * FastInverse(25), colorLower, colorUpper);
		colorv2 = Lerp(-v2.y * FastInverse(25), colorLower, colorUpper);

	}
	void ComputeNormal()
	{
		glm::vec3 e1 = v1 - v0;
		glm::vec3 e2 = v2 - v0;
		normal = glm::normalize(glm::cross(e2, e1));
	}
	vec3 Lerp(float t, vec3 a1, vec3 a2)
	{
		return a1 + t * (a2 - a1);
	}
	void Interpolate(vec3 a, vec3 b, vector<vec3>& result) {
		float t = 1.0 / result.size();
		for (int i = 0; i < result.size(); ++i) {
			result[i] = a + (t * i) * (b - a);;
		}
	}
	float GetHighestZInv(const vec3 cameraPos, const mat3 R) {
		vec3 localv0 = (v0 - cameraPos) * R;
		vec3 localv1 = (v1 - cameraPos) * R;
		vec3 localv2 = (v2 - cameraPos) * R;
		return glm::max(glm::max(1 / glm::abs(localv0[2]), 1 / glm::abs(localv1[2])), 1 / glm::abs(localv2[2]));
	}
	float FastInverse(float number)
	{
		long i;
		float x2, y;
		const float threehalfs = 1.5F;

		x2 = number * 0.5F;
		y = number;
		i = *(long*)&y;                       // evil floating point bit level hacking
		i = 0x5f3759df - (i >> 1);               // what the fuck? 
		y = *(float*)&i;
		y = y * (threehalfs - (x2 * y * y));

		return y;
	}

	vec3 fisq_normalize(vec3& v) {
		float magnitude = v.x * v.x + v.y * v.y + v.z * v.z;
		float invsqrt = FastInverse(magnitude);
		return vec3(v.x * invsqrt, v.y * invsqrt, v.z * invsqrt);
	}
	bool Render(const vec3 cameraPos, const mat3 R) {
		//has to be bigger, y-axis is upside down maxing the normals upside down
		if (glm::dot(center - cameraPos, normal) > 0)return false;

		//Calculate distance between the camera and all vertices in this polygon, return false if they are larger than drawdistance+1
		float dv0 = glm::distance(cameraPos, v0);
		if (dv0 > drawdistance + 1)return false;

		vec3 localv0 = (v0 - cameraPos) * R;
		if (localv0.z <= 0)return false;

		float dv1 = glm::distance(cameraPos, v1);
		if (dv1 > drawdistance + 1)return false;

		vec3 localv1 = (v1 - cameraPos) * R;
		if (localv1.z <= 0)return false;

		float dv2 = glm::distance(cameraPos, v2);
		if (dv2 > drawdistance + 1)return false;

		vec3 localv2 = (v2 - cameraPos) * R;
		if (localv2.z <= 0)return false;


		Pixel p0, p1, p2;
		// Fast inverse square root is faster than standard division, use that instead
		p0.zinv = FastInverse(localv0[2] * localv0[2]);
		p1.zinv = FastInverse(localv1[2] * localv1[2]);
		p2.zinv = FastInverse(localv2[2] * localv2[2]);

		//Translate to 2d screen
		p0.x = f * (localv0[0] * p0.zinv) + (SCREEN_WIDTH >> 1); //bitshift is also faster than division (focal length is screen dimensions.)
		p0.y = f * (localv0[1] * p0.zinv) + (SCREEN_HEIGHT >> 1);

		p1.x = f * (localv1[0] * p1.zinv) + (SCREEN_WIDTH >> 1);
		p1.y = f * (localv1[1] * p1.zinv) + (SCREEN_HEIGHT >> 1);


		p2.x = f * (localv2[0] * p2.zinv) + (SCREEN_WIDTH >> 1);
		p2.y = f * (localv2[1] * p2.zinv) + (SCREEN_HEIGHT >> 1);


		// Set up bounding box for polygon
		int maxX = glm::max(glm::max(p0.x, p1.x), p2.x);
		int maxY = glm::max(glm::max(p0.y, p1.y), p2.y);
		int minX = glm::min(glm::min(p0.x, p1.x), p2.x);
		int minY = glm::min(glm::min(p0.y, p1.y), p2.y);
		if (maxX < 0 || minX >= SCREEN_WIDTH || maxY < 0 || minY >= SCREEN_HEIGHT)return false;
		int boundingBoxWidth, boundingBoxHeight, x, y;
		boundingBoxWidth = maxX - minX;
		boundingBoxHeight = maxY - minY;
		float t1, t2, t3, t4;
		//Precompute values that are reused to calculate weight of vertices in a barycentric coord system
		t1 = (p1.y - p2.y);
		t2 = (p2.x - p1.x);
		t3 = (p2.y - p0.y);
		t4 = (p0.x - p2.x);
		float wDiv = (p1.y - p2.y) * (p0.x - p2.x) + (p2.x - p1.x) * (p0.y - p2.y);
		wDiv = FastInverse(wDiv * wDiv);
		float w0, w1, w2, zinv, distance;
		float draw = FastInverse((drawdistance - dimmingRange) * (drawdistance - dimmingRange));

		for (int i = 0; i < boundingBoxWidth; ++i) {
			x = minX + i;
			if (x < 0 || x >= SCREEN_WIDTH)continue;
			for (int j = 0; j < boundingBoxHeight; ++j) {
				y = minY + j;
				if (y < 0 || y >= SCREEN_HEIGHT)continue;
				//Calculate the weights of each vertex
				w0 = (t1 * (x - p2.x) + t2 * (y - p2.y)) * wDiv;
				// Check if a point in the bounding box is inside of the polygon
				if (w0 < 0)continue;
				w1 = (t3 * (x - p2.x) + t4 * (y - p2.y)) * wDiv;
				if (w1 < 0)continue;
				w2 = 1.0 - w0 - w1;
				if (w2 < 0)continue;
				else {
					zinv = w0 * p0.zinv + w1 * p1.zinv + w2 * p2.zinv;
					if (zinv > zBuffer[x][y]) {
						zBuffer[x][y] = zinv;

						distance = (w0 * dv0 + w1 * dv1 + w2 * dv2) - dimmingRange;
						if (distance < 0) {
							PutPixelSDL(screen, x, y, (w0 * colorv0 + w1 * colorv1 + w2 * colorv2) * (w0 * lightv0 + w1 * lightv1 + w2 * lightv2));
							continue;
						}
						distance = 1 - distance * draw;
						PutPixelSDL(screen, x, y, Lerp(distance, vec3(0.2823, 0.7450, 1), (w0 * colorv0 + w1 * colorv1 + w2 * colorv2) * (w0 * lightv0 + w1 * lightv1 + w2 * lightv2)));
					}
				}
			}
		}

		return true;
	}
	vector<ivec2> bres_line(int x0, int y0, int x1, int y1) {
		float dx, sx, dy, sy, err, e2;

		// Set up step directions and initial error
		dx = abs(x1 - x0);
		sx = x0 < x1 ? 1 : -1;
		dy = -abs(y1 - y0);
		sy = y0 < y1 ? 1 : -1;
		err = dx + dy;

		vector<ivec2> ret;
		// Step towards destination
		while (true) {
			ret.push_back(ivec2(x0, y0));
			if (x0 == x1 && y0 == y1) break;
			e2 = 2 * err;
			if (e2 >= dy) {
				err += dy;
				x0 += sx;
			}
			if (e2 <= dx) {
				err += dx;
				y0 += sy;
			}
		}
		return ret;
	}
	void setLights() {
		// Set light values per vertex
		vec3 temp = lightPos - v0;
		vec3 rhat = fisq_normalize(temp);

		float dot = glm::dot(rhat, normal);
		if (dot > 0) {
			vec3 directL = vec3(0, 0, 0);
			float r = glm::distance(lightPos, v0);

			directL = (lightPower * dot) * FastInverse((4 * 3.14f * r * r) * (4 * 3.14f * r * r)); // color the light
			lightv0 = directL + indirectLightPowerPerArea;
		}
		else {
			lightv0 = indirectLightPowerPerArea;
		}
		temp = lightPos - v1;
		rhat = fisq_normalize(temp);

		dot = glm::dot(rhat, normal);
		if (dot > 0) {
			vec3 directL = vec3(0, 0, 0);
			float r = glm::distance(lightPos, v1);
			directL = (lightPower * dot) * FastInverse((4 * 3.14f * r * r) * (4 * 3.14f * r * r)); // color the light
			lightv1 = directL + indirectLightPowerPerArea;
		}
		else {
			lightv1 = indirectLightPowerPerArea;
		}
		temp = lightPos - v2;
		rhat = fisq_normalize(temp);

		dot = glm::dot(rhat, normal);
		if (dot > 0) {
			vec3 directL = vec3(0, 0, 0);
			float r = glm::distance(lightPos, v2);
			directL = (lightPower * dot) * FastInverse((4 * 3.14f * r * r) * (4 * 3.14f * r * r)); // color the light
			lightv2 = directL + indirectLightPowerPerArea;
		}
		else {
			lightv2 = indirectLightPowerPerArea;
		}
	}
	void SetShadowVertices(const vector<vector<Polygon>>& polygons, vec3 lightsource) {
		if (glm::dot(lightsource - center, normal) > 0)return;
		// If the polygon faces the light source, calculate if the a polygon is shadowed by other polygons

		vec3 e1, e2, b, x;
		float t, u, v, distance;
		bool inLine = true;
		int ilow, ihigh, jlow, jhigh;
		ilow = glm::min(v0.x, lightsource.x);
		ihigh = glm::max(v0.x, lightsource.x);
		jlow = glm::min(v0.z, lightsource.z);
		jhigh = glm::max(v0.z, lightsource.z);
		vector<ivec2> line = bres_line(ilow, jlow, ihigh, jhigh);
		for (int i = 0; i < line.size(); i++) {
			e1 = polygons[line[i].x][line[i].y].v1 - polygons[line[i].x][line[i].y].v0;
			e2 = polygons[line[i].x][line[i].y].v2 - polygons[line[i].x][line[i].y].v0;
			b = v0 - polygons[line[i].x][line[i].y].v0;
			mat3 A(-(lightsource - v0), e1, e2);
			x = glm::inverse(A) * b;
			t = x[0];
			u = x[1];
			v = x[2];
			if (t > 0 && u >= 0 && v >= 0 && u + v <= 1) {
				inLine = false;
			}
		}

		if (!inLine) {
			lightv0 = indirectLightPowerPerArea;
		}

		inLine = true;
		ilow = glm::min(v1.x, lightsource.x);
		ihigh = glm::max(v1.x, lightsource.x);
		jlow = glm::min(v1.z, lightsource.z);
		jhigh = glm::max(v1.z, lightsource.z);
		line = bres_line(ilow, jlow, ihigh, jhigh);
		for (int i = 0; i < line.size(); i++) {
			e1 = polygons[line[i].x][line[i].y].v1 - polygons[line[i].x][line[i].y].v0;
			e2 = polygons[line[i].x][line[i].y].v2 - polygons[line[i].x][line[i].y].v0;
			b = v1 - polygons[line[i].x][line[i].y].v0;
			mat3 A(-(lightsource - v1), e1, e2);
			x = glm::inverse(A) * b;
			t = x[0];
			u = x[1];
			v = x[2];
			if (t > 0 && u >= 0 && v >= 0 && u + v <= 1) {
				inLine = false;
			}
		}

		if (!inLine) {
			lightv1 = indirectLightPowerPerArea;
		}
		inLine = true;
		ilow = glm::min(v2.x, lightsource.x);
		ihigh = glm::max(v2.x, lightsource.x);
		jlow = glm::min(v2.z, lightsource.z);
		jhigh = glm::max(v2.z, lightsource.z);
		line = bres_line(ilow, jlow, ihigh, jhigh);
		for (int i = 0; i < line.size(); i++) {
			e1 = polygons[line[i].x][line[i].y].v1 - polygons[line[i].x][line[i].y].v0;
			e2 = polygons[line[i].x][line[i].y].v2 - polygons[line[i].x][line[i].y].v0;
			b = v1 - polygons[line[i].x][line[i].y].v0;
			mat3 A(-(lightsource - v1), e1, e2);
			x = glm::inverse(A) * b;
			t = x[0];
			u = x[1];
			v = x[2];
			if (t > 0 && u >= 0 && v >= 0 && u + v <= 1) {
				inLine = false;
			}
		}

		if (!inLine) {
			lightv2 = indirectLightPowerPerArea;
		}

	}
};
// ----------------------------------------------------------------------------

vector<vector<float>> height(dataResX, vector<float>(dataResY));
vector<vector<float>> fertility(polyResX, vector<float>(polyResY));
vector<vector<float>> temperature(polyResX, vector<float>(polyResY));
vector<vector<Polygon>> polygons(polyResX, vector<Polygon>(polyResY));
vector<vector<bool>> renderlist(polyResX, vector<bool>(polyResY));
vector<vector<Pixel>> pixelList(SCREEN_WIDTH, vector<Pixel>(SCREEN_HEIGHT));
// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();

float FastInverse(float number);

void GenerateVectors(vector<vector<vec2>>& constants);
float Fade(float t);
float CalculatePerlin(float x, float y, vector<vector<vec2>> p);
void GenerateValues(const vector<vector<vec2>> constants, vector<vector<float>>& values, float persistance, int octave);
void GeneratePolygons(vector<vector<float>> height, vector<vector<float>> fertility, vector<vector<float>> temperature, vector<vector<Polygon>>& polygons, float sc);

void GenerateVectors(vector<vector<vec2>>& constants)
{
	for (int i = 0; i < constants.size() - 1; i++) {
		for (int j = 0; j < constants[i].size() - 1; j++) {
			constants[i][j] = vec2(float(rand() % 3 - 1), float(rand() % 3 - 1));
		}
	}
	for (int i = 0; i < constants[0].size() - 1; i++) {
		constants[constants.size() - 1][i] = constants[0][i];
	}
	for (int i = 0; i < constants.size(); i++) {
		constants[i][constants[0].size() - 1] = constants[i][0];
	}
}
float Fade(float t)
{
	return ((6 * t - 15) * t + 10) * t * t * t;
}
float Lerp(float t, float a1, float a2)
{
	return a1 + t * (a2 - a1);
}
vec3 Lerp(float t, vec3 a1, vec3 a2)
{
	return a1 + t * (a2 - a1);
}
float CalculatePerlin(float x, float y, vector<vector<vec2>> p) {
	x = fmod(x, float(p.size() - 1));
	y = fmod(y, float(p[0].size() - 1));
	float xf = (fmod(x, 1.0));
	int X = x - xf;
	float yf = (fmod(y, 1.0));
	int Y = y - yf;
	float u = Fade(xf);
	float v = Fade(yf);
	vec2 topLeft(xf, yf);
	vec2 botLeft(xf, yf - 1.0f);
	vec2 topRight(xf - 1.0f, yf);
	vec2 botRight(xf - 1.0f, yf - 1.0f);
	float floatvalue = Lerp(u, Lerp(v, glm::dot(topLeft, p[X][Y]), glm::dot(botLeft, p[X][Y + 1])), Lerp(v, glm::dot(topRight, p[X + 1][Y]), glm::dot(botRight, p[X + 1][Y + 1])));
	return (floatvalue + 1.0f) / 2.0f;
}
void GenerateValues(const vector<vector<vec2>> constants, vector<vector<float>>& values, float persistance, int octave) {
	float stepX = float(constants.size()) / float(values.size());
	float stepY = float(constants[0].size()) / float(values[0].size());
	float x, y, max, total, amplitude, frequency, lowest, highest;
	amplitude = 1;
	max = 0;
	lowest = 1.0;
	highest = 0;
	for (int a = 0; a < octave; a++) {
		max += 1.0f * amplitude;
		amplitude *= persistance;
	}
	for (int i = 0; i < values.size(); i++) {
		x = i * stepX;


		for (int j = 0; j < values[0].size(); j++) {
			y = j * stepY;
			frequency = 1.0;
			total = 0;
			amplitude = 1;
			for (int a = 0; a < octave; a++) {
				total += CalculatePerlin(x * frequency + 0.0001, y * frequency + 0.0001, constants) * amplitude;
				frequency *= 2;
				amplitude *= persistance;
			}
			values[i][j] = total / max;
			if (values[i][j] > highest) {
				highest = values[i][j];
			}
			if (values[i][j] < lowest) {
				lowest = values[i][j];
			}
		}
	}
	for (int i = 0; i < values.size(); i++) {
		for (int j = 0; j < values[0].size(); j++) {

			values[i][j] = (values[i][j] - lowest) / (highest - lowest);
		}
	}
}
void GeneratePolygons(vector<vector<float>> height, vector<vector<float>> fertility, vector<vector<float>> temperature, vector<vector<Polygon>>& polygons, float sc) {
	for (int i = 0; i < height.size() - 1; i++) {
		for (int j = 0; j < height[0].size() - 1; j++) {
			polygons[i * 2][j] = Polygon(vec3(i * sc, 5 * height[i][j] - 1, j * sc), vec3((i + 1) * sc, 5 * height[i + 1][j + 1] - 1, (j + 1) * sc), vec3((i + 1) * sc, 5 * height[i + 1][j] - 1, j * sc), fertility[i * 2][j], temperature[i * 2][j]);
			polygons[i * 2 + 1][j] = Polygon(vec3(i * sc, 5 * height[i][j] - 1, j * sc), vec3(i * sc, 5 * height[i][j + 1] - 1, (j + 1) * sc), vec3((i + 1) * sc, 5 * height[i + 1][j + 1] - 1, (j + 1) * sc), fertility[i * 2 + 1][j], temperature[i * 2 + 1][j]);
		}
	}
}

//WITH AN APPEARANCE OF HIT ARTIST AND MUSIC PRODUCER, QUAKE III?
float FastInverse(float number)
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;

	x2 = number * 0.5F;
	y = number;
	i = *(long*)&y;                       // evil floating point bit level hacking
	i = 0x5f3759df - (i >> 1);               // what the fuck? 
	y = *(float*)&i;
	y = y * (threehalfs - (x2 * y * y));

	return y;
}

int main(int argc, char* argv[])
{
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.
	float* depthBuffer = (float*)malloc(sizeof(float) * SCREEN_HEIGHT * SCREEN_WIDTH);
	vector<vector<vec2>> h(5, vector<vec2>(5));
	vector<vector<vec2>> f(4, vector<vec2>(4));
	vector<vector<vec2>> t(4, vector<vec2>(4));
	GenerateVectors(h);
	GenerateVectors(f);
	GenerateVectors(t);
	GenerateValues(h, height, h_persistence, h_octave);
	GenerateValues(f, fertility, f_persistence, f_octave);
	GenerateValues(t, temperature, t_persistence, t_octave);
	GeneratePolygons(height, fertility, temperature, polygons, scale);
	for (int i = 0; i < polygons.size(); ++i) {
		for (int j = 0; j < polygons[i].size(); ++j) {
			polygons[i][j].SetShadowVertices(polygons, lightPos);
		}
	}

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;
	//cout << "fps: " << 1000 / dt << endl;
	Uint8* keystate = SDL_GetKeyState(0);

	if (keystate[SDLK_UP])
		;

	if (keystate[SDLK_DOWN])
		;

	if (keystate[SDLK_RIGHT])
		;

	if (keystate[SDLK_LEFT])
		;

	if (keystate[SDLK_RSHIFT]) {
		//reset camera
		cameraPos = vec3(0, 0, -3.001);
		R = mat3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
		yaw = 0;
	}

	if (keystate[SDLK_RCTRL]) {
		// change light position to current cameraPos
		lightPos = cameraPos;
		// Light can't be out of bounds of the landscape
		if (lightPos[0] < 0)
			lightPos[0] = 0;
		if (lightPos[0] >= dataResX * scale)
			lightPos[0] = dataResX * scale - 1;
		if (lightPos[2] < 0)
			lightPos[2] = 0;
		if (lightPos[2] >= dataResY * scale)
			lightPos[2] = dataResY * scale - 1;
		cout << "Relighting...";
		for (int i = 0; i < polygons.size(); ++i) {
			for (int j = 0; j < polygons[i].size(); ++j) {
				polygons[i][j].setLights();
				polygons[i][j].SetShadowVertices(polygons, lightPos);
			}
		}
		cout << " done" << endl;
	}

	if (keystate[SDLK_w]) {
		// Move camera forward
		cameraPos[0] += 0.1 * (-sin(yaw * turnRate));
		cameraPos[2] += 0.1 * cos(yaw * turnRate);
	}

	else if (keystate[SDLK_s]) {
		// Move camera backward
		cameraPos[0] -= 0.1 * (-sin(yaw * turnRate));
		cameraPos[2] -= 0.1 * cos(yaw * turnRate);
	}

	if (keystate[SDLK_d]) {
		// Move camera to the right
		cameraPos[0] += 0.1 * cos(yaw * turnRate);
		cameraPos[2] += 0.1 * sin(yaw * turnRate);
	}

	else if (keystate[SDLK_a]) {
		// Move camera to the left
		cameraPos[0] -= 0.1 * cos(yaw * turnRate);
		cameraPos[2] -= 0.1 * sin(yaw * turnRate);
	}

	if (keystate[SDLK_e]) {
		// Rotate camera to the right
		yaw -= 3.14 / 36;
		R[0][0] = cos(yaw * turnRate);
		R[0][2] = sin(yaw * turnRate);
		R[2][0] = -sin(yaw * turnRate);
		R[2][2] = cos(yaw * turnRate);
	}

	else if (keystate[SDLK_q]) {
		// Rotate camera to the left
		yaw += 3.14 / 36;
		R[0][0] = cos(yaw * turnRate);
		R[0][2] = sin(yaw * turnRate);
		R[2][0] = -sin(yaw * turnRate);
		R[2][2] = cos(yaw * turnRate);
	}
	if (keystate[SDLK_SPACE]) {
		// Rotate camera to the left
		cameraPos[1] -= 0.1;
	}
	else if (keystate[SDLK_LCTRL]) {
		// Rotate camera to the left
		cameraPos[1] += 0.1;
	}
}
void Draw()
{

	SDL_FillRect(screen, 0, 0x47beff);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	// Refill depth buffer
	for (int i = 0; i < SCREEN_WIDTH; i++)
		for (int j = 0; j < SCREEN_HEIGHT; j++)
			zBuffer[i][j] = 0;

	int culled = 0;

	// Calls the draw function for each polygon
	for (int i = 0; i < polygons.size(); ++i) {
		for (int j = 0; j < polygons[i].size(); ++j) {
			if (!polygons[i][j].Render(cameraPos, R))
				culled++;
		}
	}

	//cout << float(culled) * FastInverse(polygons.size() * polygons[0].size() * polygons.size() * polygons[0].size()) << "% culled " << endl;
	//cout << "pos: [" << cameraPos[0] << ", " << cameraPos[1] << ", " << cameraPos[2] << "]" << endl;
	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}