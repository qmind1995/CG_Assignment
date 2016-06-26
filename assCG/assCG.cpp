// Tran Quoc Tri
//	51304367

#include "stdafx.h"
#include <iostream>
#include <windows.h>
#include <gl.h>
#include <glut.h>
#include <math.h>
#include <stdio.h>										// Standard I/O header 
#include <stdlib.h>

typedef	struct
{
	GLubyte	* imageData;									// Image Data (Up To 32 Bits)
	GLuint	bpp;											// Image Color Depth In Bits Per Pixel
	GLuint	width;											// Image Width
	GLuint	height;											// Image Height
	GLuint	texID;											// Texture ID Used To Select A Texture
	GLuint	type;											// Image Type (GL_RGB, GL_RGBA)
} Texture;

typedef struct
{
	GLubyte Header[12];									// TGA File Header
} TGAHeader;

typedef struct
{
	GLubyte		header[6];								// First 6 Useful Bytes From The Header
	GLuint		bytesPerPixel;							// Holds Number Of Bytes Per Pixel Used In The TGA File
	GLuint		imageSize;								// Used To Store The Image Size When Setting Aside Ram
	GLuint		temp;									// Temporary Variable
	GLuint		type;
	GLuint		Height;									//Height of Image
	GLuint		Width;									//Width ofImage
	GLuint		Bpp;									// Bits Per Pixel
} TGA;


TGAHeader tgaheader;									// TGA header
TGA tga;										// TGA image data

Texture   floorTex;

GLubyte uTGAcompare[12] = { 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };	// Uncompressed TGA Header


bool LoadTGA(Texture * texture, char * filename)				// Load a TGA file
{
	FILE * fTGA;												// File pointer to texture file
	fTGA = fopen(filename, "rb");								// Open file for reading

	if (fTGA == NULL)											// If it didn't open....
	{
		return false;														// Exit function
	}

	if (fread(&tgaheader, sizeof(TGAHeader), 1, fTGA) == 0)					// Attempt to read 12 byte header from file
	{
		if (fTGA != NULL)													// Check to seeiffile is still open
		{
			fclose(fTGA);													// If it is, close it
		}
		return false;														// Exit function
	}

	// an Uncompressed TGA image
	if (fread(tga.header, sizeof(tga.header), 1, fTGA) == 0)					// Read TGA header
	{
		if (fTGA != NULL)													// if file is still open
		{
			fclose(fTGA);													// Close it
		}
		return false;														// Return failular
	}

	texture->width = tga.header[1] * 256 + tga.header[0];					// Determine The TGA Width	(highbyte*256+lowbyte)
	texture->height = tga.header[3] * 256 + tga.header[2];					// Determine The TGA Height	(highbyte*256+lowbyte)
	texture->bpp = tga.header[4];										// Determine the bits per pixel
	tga.Width = texture->width;										// Copy width into local structure						
	tga.Height = texture->height;										// Copy height into local structure
	tga.Bpp = texture->bpp;											// Copy BPP into local structure

	if ((texture->width <= 0) || (texture->height <= 0) || ((texture->bpp != 24) && (texture->bpp != 32)))	// Make sure all information is valid
	{
		if (fTGA != NULL)													// Check if file is still open
		{
			fclose(fTGA);													// If so, close it
		}
		return false;														// Return failed
	}

	if (texture->bpp == 24)													// If the BPP of the image is 24...
		texture->type = GL_RGB;											// Set Image type to GL_RGB
	else																	// Else if its 32 BPP
		texture->type = GL_RGBA;											// Set image type to GL_RGBA

	tga.bytesPerPixel = (tga.Bpp / 8);									// Compute the number of BYTES per pixel
	tga.imageSize = (tga.bytesPerPixel * tga.Width * tga.Height);		// Compute the total amout ofmemory needed to store data
	texture->imageData = (GLubyte *)malloc(tga.imageSize);					// Allocate that much memory

	if (texture->imageData == NULL)											// If no space was allocated
	{
		fclose(fTGA);														// Close the file
		return false;														// Return failed
	}

	if (fread(texture->imageData, 1, tga.imageSize, fTGA) != tga.imageSize)	// Attempt to read image data
	{
		if (texture->imageData != NULL)										// If imagedata has data in it
		{
			free(texture->imageData);										// Delete data from memory
		}
		fclose(fTGA);														// Close file
		return false;														// Return failed
	}

	// switch R and B
	for (int i = 0; i < tga.imageSize; i += tga.bytesPerPixel)
	{
		GLubyte temp = texture->imageData[i];
		texture->imageData[i] = texture->imageData[i + 2];
		texture->imageData[i + 2] = temp;
	}


	fclose(fTGA);															// Close file
	return true;															// All went well, continue on
}
bool on_off=true;
int		screenWidth = 600;
int		screenHeight = 600;
#define	COLORNUM		17//14
#define PI			3.1415926

bool	bWireFrame = false;
float	baseRadius = 0.6;
float	baseHeight = 0.2;
float	baseRotateStep = 5;

float	camRotateStep = 5;
float	camAngle = 135;
float	camUpAngle = 0;

float	cylRadius = 0.3;
float	cylHeight = 1.0;
float   cylMaxScaleY = 2.0;
float	cylScaleStep = 0.05;

bool anima=false;
bool trend = true; //true = incre
bool fourviewmode = false;
bool plustrend = true; // true trend

float	body1SizeX = 3;
float	body1SizeY = 0.2;
float	body1SizeZ = 0.8;
int selection = 0;
bool isR = false;
bool isG = false;
bool isB = false;
using namespace std;

float	ColorArr[COLORNUM][3] = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 },
{ 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 },
{ 0.3, 0.3, 0.3 }, { 0.5, 0.5, 0.5 }, { 0.9, 0.9, 0.9 },
{ 1.0, 0.5, 0.5 }, { 0.5, 1.0, 0.5 }, { 0.5, 0.5, 1.0 },
{ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, { 8.0 / 255, 39.0 / 255, 99.0 / 255 }, {57.0/255,0.0,201.0/255} };

class Point3
{
public:
	float x, y, z;
	void set(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}
	void set(Point3& p)
	{
		x = p.x; y = p.y; z = p.z;
	}
	Point3() { x = y = z = 0; }
	Point3(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}

};

class VertexID
{
public:
	int		vertIndex;
	int		colorIndex;
};

class Vector3{
public:
	float	x, y, z;
	void set(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}
	void set(Vector3& v)
	{
		x = v.x; y = v.y; z = v.z;
	}
	void flip()
	{
		x = -x; y = -y; z = -z;
	}
	//void setDiff(Point3& a, Point3&
	void normalize();
	Vector3() { x = y = z = 0; }
	Vector3(float dx, float dy, float dz)
	{
		x = dx; y = dy; z = dz;
	}
	Vector3(Vector3& v)
	{
		x = v.x; y = v.y; z = v.z;
	}
	Vector3 cross(Vector3 b);
	float dot(Vector3 b);
};

class Face
{
public:
	int		nVerts;
	Vector3 facenorm;
	VertexID*	vert;

	Face()
	{
		nVerts = 0;
		vert = NULL;
	}
	~Face()
	{
		if (vert != NULL)
		{
			delete[] vert;
			vert = NULL;
		}
		nVerts = 0;
	}
};

class Mesh
{
public:
	int		numVerts;
	Point3*		pt;

	float	slideX, slideY, slideZ;
	float	rotateX, rotateY, rotateZ;
	float	scaleX, scaleY, scaleZ;

	int		numFaces;
	Face*		face;
public:
	Mesh()
	{
		numVerts = 0;
		pt = NULL;
		numFaces = 0;
		face = NULL;
	}
	~Mesh()
	{
		if (pt != NULL)
		{
			delete[] pt;
		}
		if (face != NULL)
		{
			delete[] face;
		}
		numVerts = 0;
		numFaces = 0;
	}
	void DrawWireframe();
	void DrawColor();
	void CreateCuboid(float	fSizeX, float fSizeY, float	fSizeZ);
	void CreateTetrahedron();
	void CreateCylinder(int nSegment, float fHeight, float fRadius);
	void Creat4(int nSegment, float r1, float r2, float fLength, float fHeight);
	void CreateSlipper(float r1, float r2, float len,float height);
	void CreateSphere(int nSlice, int nStack, float radius);
	void CreateCube(float	fSize);
	void CreateFrame(float len, float width, float height);
	void SetColor(int colorIdx);
	void Paving();
	void CalculateFacesNorm();
	void Draw();
};


static void myShadowMatrix(float ground[4], float light[4])
{
	float  dot;
	float  shadowMat[4][4];

	dot = ground[0] * light[0] +
		ground[1] * light[1] +
		ground[2] * light[2] +
		ground[3] * light[3];

	shadowMat[0][0] = dot - light[0] * ground[0];
	shadowMat[1][0] = 0.0 - light[0] * ground[1];
	shadowMat[2][0] = 0.0 - light[0] * ground[2];
	shadowMat[3][0] = 0.0 - light[0] * ground[3];

	shadowMat[0][1] = 0.0 - light[1] * ground[0];
	shadowMat[1][1] = dot - light[1] * ground[1];
	shadowMat[2][1] = 0.0 - light[1] * ground[2];
	shadowMat[3][1] = 0.0 - light[1] * ground[3];

	shadowMat[0][2] = 0.0 - light[2] * ground[0];
	shadowMat[1][2] = 0.0 - light[2] * ground[1];
	shadowMat[2][2] = dot - light[2] * ground[2];
	shadowMat[3][2] = 0.0 - light[2] * ground[3];

	shadowMat[0][3] = 0.0 - light[3] * ground[0];
	shadowMat[1][3] = 0.0 - light[3] * ground[1];
	shadowMat[2][3] = 0.0 - light[3] * ground[2];
	shadowMat[3][3] = dot - light[3] * ground[3];

	glMultMatrixf((const GLfloat*)shadowMat);
}

class Camera
{
private:
	Point3       eye,lookat;
	Vector3      u, v, n;
	double       viewAngle, aspect, nearDist, farDist;
	void         setModelViewMatrix();
public:
	Camera(){};
	void  set(Point3 Eye, Point3 look, Vector3 up);
	void  flyCam(float);
	void  slideXZ(float delU, float delV, float delN);
	void  setShape(float vAng, float asp, float nearD, float farD);
	void rotaY();
};

Vector3 Vector3::cross(Vector3 b){
	Vector3 c(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);
	return c;
}

void Camera::setModelViewMatrix()
{
	float         m[16];
	Vector3   eVec(eye.x, eye.y, eye.z);
	m[0] = u.x; m[4] = u.y; m[8] = u.z; m[12] = -eVec.dot(u);
	m[1] = v.x; m[5] = v.y; m[9] = v.z; m[13] = -eVec.dot(v);
	m[2] = n.x; m[6] = n.y; m[10] = n.z; m[14] = -eVec.dot(n);
	m[3] = 0; m[7] = 0; m[11] = 0; m[15] = 1.0;
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(m);
}

void Camera::set(Point3 Eye, Point3 look, Vector3 up){
	lookat.set(look);
	eye.set(Eye);
	n.set(eye.x - look.x, eye.y - look.y, eye.z - look.z);
	u.set(up.cross(n));
	n.normalize();
	u.normalize();
	v.set(n.cross(u));
	setModelViewMatrix();
}

void Camera::slideXZ(float delU, float delV, float delN)
{
	eye.x += delU*u.x + delV*v.x + delN*n.x;
	eye.z += delU*u.z + delV*v.z + delN*n.z;
	
	set(eye, lookat, Vector3(0, 1, 0));
	setModelViewMatrix();
}

void Camera::rotaY(){
	float d = sqrt(pow(eye.x - lookat.x, 2)  + pow(eye.z - lookat.z, 2));
	eye.x = d*cos(camAngle*PI / 180);
	eye.z = d*sin(camAngle*PI / 180);
	set(eye, lookat, Vector3(0, 1, 0));
	setModelViewMatrix();
}

void Camera::flyCam(float step){
	eye.y += step;
	set(eye, lookat, Vector3(0, 1, 0));
	setModelViewMatrix();
}

void  Camera::setShape(float vAng, float asp, float nearD, float farD)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	viewAngle = vAng;
	aspect = asp;
	nearDist = nearD;
	farDist = farD;
	gluPerspective(viewAngle, aspect, nearDist, farDist);
	setModelViewMatrix();
}
Camera cam;

float Vector3::dot(Vector3 b){
	return x*b.x + y*b.y + z*b.z;
}

void Vector3::normalize(){
	float temp = sqrt(x*x + y*y + z*z);
	x = x / temp;
	y = y / temp;
	z = z / temp;
}

void Mesh::DrawColor()
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (int f = 0; f < numFaces; f++)
	{
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++)
		{
			int		iv = face[f].vert[v].vertIndex;
			int		ic = face[f].vert[v].colorIndex;

			glColor3f(ColorArr[ic][0], ColorArr[ic][1], ColorArr[ic][2]);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

void Mesh::SetColor(int colorIdx){
	for (int f = 0; f < numFaces; f++)
	{
		for (int v = 0; v < face[f].nVerts; v++)
		{
			face[f].vert[v].colorIndex = colorIdx;
		}
	}
}

void Mesh::DrawWireframe()
{

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int f = 0; f < numFaces; f++)
	{
		
		glBegin(GL_POLYGON);
		
		for (int v = 0; v < face[f].nVerts; v++)
		{
			//glLineWidth(5.0);
			int		iv = face[f].vert[v].vertIndex;

			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Mesh::CreateCylinder(int nSegment, float fHeight, float fRadius){
	int i=0, j=0;
	numVerts = nSegment + 1;
	pt = new Point3[numVerts * 2];
	for (i = 0; i < nSegment; i++)
		pt[i].set(fRadius*sin(2 * i * PI / nSegment), 0, fRadius*cos(2 * i * PI / nSegment));

	for (i = nSegment; i < 2 * nSegment; i++)
		pt[i].set(fRadius*sin(2 * i * PI / nSegment), fHeight, fRadius*cos(2 * i * PI / nSegment));

	pt[2 * nSegment].set(0, 0, 0);
	pt[2 * nSegment + 1].set(0, fHeight, 0);

	numFaces = nSegment * 3;
	face = new Face[numFaces];
	//ve da giac cuoi cung
	face[nSegment - 1].nVerts = 3;
	face[nSegment - 1].vert = new VertexID[face[nSegment - 1].nVerts];
	face[nSegment - 1].vert[0].vertIndex = nSegment - 1;
	face[nSegment - 1].vert[1].vertIndex = 0;
	face[nSegment - 1].vert[2].vertIndex = 2 * nSegment;

	face[2 * nSegment - 1].nVerts = 3;
	face[2 * nSegment - 1].vert = new VertexID[face[nSegment - 1].nVerts];
	face[2 * nSegment - 1].vert[0].vertIndex = 2 * nSegment - 1;
	face[2 * nSegment - 1].vert[1].vertIndex = nSegment;
	face[2 * nSegment - 1].vert[2].vertIndex = 2 * nSegment + 1;

	face[3 * nSegment - 1].nVerts = 4;
	face[3 * nSegment - 1].vert = new VertexID[face[3 * nSegment - 1].nVerts];
	face[3 * nSegment - 1].vert[0].vertIndex = nSegment - 1;
	face[3 * nSegment - 1].vert[1].vertIndex = 0;
	face[3 * nSegment - 1].vert[2].vertIndex = nSegment;
	face[3 * nSegment - 1].vert[3].vertIndex = 2 * nSegment - 1;

	for (i = 0; i < nSegment - 1; i++)
	{
		//mat tron 1
		face[i].nVerts = 3;
		face[i].vert = new VertexID[face[i].nVerts];
		face[i].vert[0].vertIndex = i;
		face[i].vert[1].vertIndex = i + 1;
		face[i].vert[2].vertIndex = 2 * nSegment;// tam duong tron 1

		//mat tron 2
		face[i + nSegment].nVerts = 3;
		face[i + nSegment].vert = new VertexID[face[i + nSegment].nVerts];

		face[i + nSegment].vert[0].vertIndex = i + nSegment;
		face[i + nSegment].vert[1].vertIndex = i + nSegment + 1;
		face[i + nSegment].vert[2].vertIndex = 2 * nSegment + 1;// tam duong tron 2

		//mat ben
		face[i + 2 * nSegment].nVerts = 4;
		face[i + 2 * nSegment].vert = new VertexID[face[i + 2 * nSegment].nVerts];

		face[i + 2 * nSegment].vert[0].vertIndex = i;
		face[i + 2 * nSegment].vert[1].vertIndex = i + 1;
		face[i + 2 * nSegment].vert[2].vertIndex = i + nSegment + 1;
		face[i + 2 * nSegment].vert[3].vertIndex = i + nSegment;

	}
}

void Mesh::CreateCuboid(float fSizeX, float fSizeY, float fSizeZ){
	int i=0;
	numVerts = 8;
	pt = new Point3[numVerts];
	pt[0].set(-fSizeX, fSizeY, fSizeZ);
	pt[1].set(fSizeX, fSizeY, fSizeZ);
	pt[2].set(fSizeX, fSizeY, -fSizeZ);
	pt[3].set(-fSizeX, fSizeY, -fSizeZ);
	pt[4].set(-fSizeX, -fSizeY, fSizeZ);
	pt[5].set(fSizeX, -fSizeY, fSizeZ);
	pt[6].set(fSizeX, -fSizeY, -fSizeZ);
	pt[7].set(-fSizeX, -fSizeY, -fSizeZ);

	numFaces = 6;
	face = new Face[numFaces];
	//Left face
	face[0].nVerts = 4;
	face[0].vert = new VertexID[face[0].nVerts];
	face[0].vert[0].vertIndex = 1;
	face[0].vert[1].vertIndex = 5;
	face[0].vert[2].vertIndex = 6;
	face[0].vert[3].vertIndex = 2;

	//Right face
	face[1].nVerts = 4;
	face[1].vert = new VertexID[face[1].nVerts];
	face[1].vert[0].vertIndex = 0;
	face[1].vert[1].vertIndex = 3;
	face[1].vert[2].vertIndex = 7;
	face[1].vert[3].vertIndex = 4;

	//top face
	face[2].nVerts = 4;
	face[2].vert = new VertexID[face[2].nVerts];
	face[2].vert[0].vertIndex = 0;
	face[2].vert[1].vertIndex = 1;
	face[2].vert[2].vertIndex = 2;
	face[2].vert[3].vertIndex = 3;

	//bottom face
	face[3].nVerts = 4;
	face[3].vert = new VertexID[face[3].nVerts];
	face[3].vert[0].vertIndex = 7;
	face[3].vert[1].vertIndex = 6;
	face[3].vert[2].vertIndex = 5;
	face[3].vert[3].vertIndex = 4;

	//near face
	face[4].nVerts = 4;
	face[4].vert = new VertexID[face[4].nVerts];
	face[4].vert[0].vertIndex = 4;
	face[4].vert[1].vertIndex = 5;
	face[4].vert[2].vertIndex = 1;
	face[4].vert[3].vertIndex = 0;

	//Far face
	face[5].nVerts = 4;
	face[5].vert = new VertexID[face[5].nVerts];
	face[5].vert[0].vertIndex = 3;
	face[5].vert[1].vertIndex = 2;
	face[5].vert[2].vertIndex = 6;
	face[5].vert[3].vertIndex = 7;
}

void Mesh::CreateSlipper(float r1, float r2, float len,float height){
	int nslice = 10;
	numVerts = (nslice+1) * 4;
	float delta = 180.0f / nslice;
	pt = new Point3[numVerts];
	for (int i = 0; i <= nslice; i++){
		float x1 = len / 2 + r1* sin(delta*i*PI / 180);
		float z1 = r1*cos(delta*i*PI / 180);
		float x2 = -(len / 2 + r2*sin(delta*i*PI / 180));
		float z2 = -r2*cos(delta*i*PI / 180);
		pt[i].set(x1, height / 2, z1); 
		pt[i + nslice + 1].set(x2, height / 2, z2);
		pt[i + 2 * nslice + 2].set(x1, -height / 2, z1);
		pt[i + 3 * nslice + 3].set(x2, -height / 2, z2);
	}
	numFaces = 2 * nslice +4;
	face = new Face[numFaces];

	for (int i = 0; i <2 * nslice+1; i++){
		face[i].nVerts = 4;
		face[i].vert = new VertexID[face[i].nVerts];
		face[i].vert[0].vertIndex = i;
		face[i].vert[3].vertIndex = i + 1;
		face[i].vert[2].vertIndex = i + 2 * nslice + 3;
		face[i].vert[1].vertIndex = i + 2 * nslice + 2;
	}

	face[2 * nslice + 1].nVerts = 4;
	face[2 * nslice + 1].vert = new VertexID[face[2 * nslice + 1].nVerts];
	face[2 * nslice + 1].vert[0].vertIndex = 0;
	face[2 * nslice + 1].vert[1].vertIndex = 2*nslice+1;
	face[2 * nslice + 1].vert[2].vertIndex = 4*nslice+3;
	face[2 * nslice + 1].vert[3].vertIndex = 2*nslice+2;

	face[2 * nslice+2].nVerts = (nslice + 1) * 2;
	face[2 * nslice+2].vert = new VertexID[face[2 * nslice+2].nVerts];
	for (int i = 0; i < 2 * nslice + 2; i++){
		face[2 * nslice+2].vert[i].vertIndex = i;
	}

	face[2 * nslice + 3].nVerts = (nslice + 1) * 2;
	face[2 * nslice+3].vert = new VertexID[face[2 * nslice+3].nVerts];
	for (int i = 0; i < 2 * nslice + 2; i++){
		face[2 * nslice+3].vert[i].vertIndex = i+2*nslice+2;
	}

}

void Mesh::Paving(){
	
	int i=0;
	float x = sin(60 * PI / 180);
	numVerts = 7;

	pt = new Point3[numVerts];
	pt[0].set(0, 0, 0);
	pt[1].set(0, 0, 1);
	pt[2].set(x, 0, 1.5);
	pt[3].set(x, 0, 0.5);
	
	pt[4].set(2*x, 0, 1);
	pt[5].set(2*x, 0, 0);
	pt[6].set(x, 0, -0.5);

	numFaces = 3;
	face = new Face[numFaces];

	face[0].nVerts = 4;
	face[0].vert = new VertexID[face[0].nVerts];
	face[0].vert[0].vertIndex = 0;
	face[0].vert[1].vertIndex = 1;
	face[0].vert[2].vertIndex = 2;
	face[0].vert[3].vertIndex = 3;
	for (i = 0; i<face[0].nVerts; i++)
		face[0].vert[i].colorIndex = 15;

	face[1].nVerts = 4;
	face[1].vert = new VertexID[face[1].nVerts];
	face[1].vert[0].vertIndex = 2;
	face[1].vert[1].vertIndex = 4;
	face[1].vert[2].vertIndex = 5;
	face[1].vert[3].vertIndex = 3;
	for (i = 0; i<face[1].nVerts; i++)
		face[1].vert[i].colorIndex = 2;

	face[2].nVerts = 4;
	face[2].vert = new VertexID[face[2].nVerts];
	face[2].vert[0].vertIndex = 0;
	face[2].vert[1].vertIndex = 3;
	face[2].vert[2].vertIndex = 5;
	face[2].vert[3].vertIndex = 6;
	for (i = 0; i < face[2].nVerts; i++)
		face[2].vert[i].colorIndex = 14;
}

void Mesh::CalculateFacesNorm(){
	for (int f = 0; f < numFaces; f++){
		Point3  p1 =  pt[face[f].vert[0].vertIndex];
		Point3  p2 = pt[face[f].vert[1].vertIndex];
		Point3  p3 = pt[face[f].vert[2].vertIndex];

		Vector3 v1 = Vector3(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
		Vector3 v2 = Vector3(p1.x - p3.x, p1.y - p3.y, p1.z - p3.z);
		face[f].facenorm = v1.cross(v2);
		face[f].facenorm.normalize();
	}
}

void Mesh::Draw(){
	for (int f = 0; f < numFaces; f++){
		glBegin(GL_POLYGON);
		for (int v = 0; v < face[f].nVerts; v++){
			int		iv = face[f].vert[v].vertIndex;
			glNormal3f(face[f].facenorm.x, face[f].facenorm.y, face[f].facenorm.z);
			glVertex3f(pt[iv].x, pt[iv].y, pt[iv].z);
		}
		glEnd();
	}
}

Mesh	base;
Mesh	cyl;
Mesh	body1;
Mesh	pav;
Mesh	cottar1, cottar2, cottar3;
Mesh	arm1,arm2,arm3;
Mesh	slide;
Mesh	frame1, frame2, frame3;
GLfloat	diffuse[] = { 0.5f, 0.5f, 0.5f, 1.0f };
GLfloat	specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat	ambient[] = { 0, 1, 0, 1.0f };

GLfloat	shiness = 100;

GLfloat arm1ambient[] = { 1, 1, 0, 0.5f };

GLfloat arm2ambient[] = { 0.0, 0.0, 1.0, 1.0f };

GLfloat slideambient[] = { 0.25, 0.25, 0.25, 0.5f };

GLfloat reddiff[] = { 1, 0, 0, 0 };
GLfloat greendiff[] = { 0.0, 1, 0, 0};
GLfloat bluediff[] = { 0, 0, 1, 0 };
GLfloat pdiff[] = { 0.5, 0, 0.5, 0 };

void loadTextures(void)	{
	bool status = LoadTGA(&floorTex, "marble.tga");
	if (status) {
		glGenTextures(1, &floorTex.texID);
		glBindTexture(GL_TEXTURE_2D, floorTex.texID);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, floorTex.width,
			floorTex.height, 0,
			GL_RGB, GL_UNSIGNED_BYTE, floorTex.imageData);

		if (floorTex.imageData)
			free(floorTex.imageData);
	}
}

void processTimer(int value){
	if (anima){
		
		if (cyl.scaleY > cylMaxScaleY || cyl.scaleY < 1) trend = !trend;
		if (trend) cyl.scaleY += cylScaleStep;
		else cyl.scaleY -= cylScaleStep;
		
		cottar1.rotateY += baseRotateStep;
		if (cottar1.rotateY>360)
			cottar1.rotateY -= 360;

		base.rotateY += baseRotateStep;
		if (base.rotateY > 360)
			base.rotateY -= 360;


		glutTimerFunc(50, processTimer, value);
		glutPostRedisplay();
	}
}

void keyboardFunc(int key, int x, int y){
	switch (key){
	case GLUT_KEY_RIGHT:
		camAngle -= 5;
		if (camAngle<360)
			camAngle += 360;
		cam.rotaY();
		break;
	case GLUT_KEY_LEFT:
		camAngle += 5;
		if (camAngle>360)
			camAngle -= 360;
		cam.rotaY();
		break;
	case GLUT_KEY_UP:

		cam.flyCam(1);
		break;
	case GLUT_KEY_DOWN:
		cam.flyCam(-1);
		break;
	}
	glutPostRedisplay();
}

void myKeyboard(unsigned char key, int x, int y){
	float	fRInc=0;
	float	fAngle=0;
	switch (key)
	{
	case '1':
		base.rotateY += baseRotateStep;
		if (base.rotateY > 360)
			base.rotateY -= 360;
		break;
	case '2':
		base.rotateY -= baseRotateStep;
		if (base.rotateY < 0)
			base.rotateY += 360;
		break;
	case '3':
		cottar1.rotateY += baseRotateStep;
		if (cottar1.rotateY>360)
			cottar1.rotateY -= 360;
		break;
	case '4':
		cottar1.rotateY -= baseRotateStep;
		if (cottar1.rotateY<360)
			cottar1.rotateY += 360;
		break;
	case 'w':
	case 'W':
		bWireFrame = !bWireFrame;
		break;
	case 'l':
	case'L':
		cyl.scaleY += cylScaleStep;
		break;
	case'x':
	case'X':
		cyl.scaleY -= cylScaleStep;
		break;
	case'-':
		cam.slideXZ(0, 0, -0.2);
		break;
	case'+':
		cam.slideXZ(0, 0, 0.2);
		break;
	case'v':
	case'V':
		fourviewmode = !fourviewmode;
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glOrtho(-4, 4, -4, 4, -1000, 1000);
		gluLookAt(0.0, 20.0, 0.0, 0, 0, 0, 0, 0 , 1);
		break;
	case'a':
	case'A':
		anima = !anima;
		if (anima) glutTimerFunc(50, processTimer, 1);
		break;
	case 'd':
	case'D':
		on_off = !on_off;
		break;
	case 'r':
	case 'R':
		if (selection == 1){
			for (int i = 0; i < 4; i++)
				arm1ambient[i] = reddiff[i];
		}
		else if (selection == 2){
			for (int i = 0; i < 4; i++)
				arm2ambient[i] = reddiff[i];
		}
		else if (selection == 3){
			for (int i = 0; i < 4; i++)
				slideambient[i] = reddiff[i];
		}
		break;
	case 'g':
	case 'G':
		if (selection == 1){
			for (int i = 0; i < 4; i++)
				arm1ambient[i] = greendiff[i];
		}
		else if (selection == 2){
			for (int i = 0; i < 4; i++)
				arm2ambient[i] = greendiff[i];
		}
		else if (selection == 3){
			for (int i = 0; i < 4; i++)
				slideambient[i] = greendiff[i];
		}
		break;
	case 'b':
	case 'B':
		if (selection == 1){
			for (int i = 0; i < 4; i++)
				arm1ambient[i] = bluediff[i];
		}
		else if (selection == 2){
			for (int i = 0; i < 4; i++)
				arm2ambient[i] = bluediff[i];
		}
		else if (selection == 3){
			for (int i = 0; i < 4; i++)
				slideambient[i] = bluediff[i];
		}
		break;
	
	case'P':
	case'p':
		if (selection == 1){
			for (int i = 0; i < 4; i++)
				arm1ambient[i] = pdiff[i];
		}
		else if (selection == 2){
			for (int i = 0; i < 4; i++)
				arm2ambient[i] = pdiff[i];
		}
		else if (selection == 3){
			for (int i = 0; i < 4; i++)
				slideambient[i] = pdiff[i];
		}
		break;
	}
	glutPostRedisplay();
}

void onMouse(int button, int state, int x, int y) {
	if (state != GLUT_DOWN)
		return;
	GLbyte color[4];
	GLfloat depth;
	GLuint index;

	glReadPixels(x, 600 - y - 1, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, color);
	glReadPixels(x, 600 - y - 1, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
	glReadPixels(x, 600 - y - 1, 1, 1, GL_STENCIL_INDEX, GL_UNSIGNED_INT, &index);
	selection = index;
	/*printf("Clicked on pixel %d, %d, color %02hhx%02hhx%02hhx%02hhx, depth %f, stencil index %u\n",
		x, y, color[0], color[1], color[2], color[3], depth, index);*/
	glutPostRedisplay();
}

void setupMaterial(float ambient[], float diffuse[], float specular[], float shiness)
{
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shiness);
}

void drawBase(){
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);

	if (bWireFrame)
		base.DrawWireframe();
	else{
		GLfloat baseambient[] = { 1.0, 0.0, 0.0, 1.0f };
		setupMaterial(baseambient, diffuse, specular, shiness);
		base.Draw();
	}
	glPopMatrix();
}

void drawBody1(){
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(0, cylHeight*cyl.scaleY + baseHeight*2, 0);
	body1.SetColor(7);
	if (bWireFrame)
		body1.DrawWireframe();
	else{

		setupMaterial(ambient, diffuse, specular, shiness);
		body1.Draw();
	}
		

	glPopMatrix();
}

void drawBody2(){
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	
	glTranslated(-0.25 * body1SizeX, body1SizeY*2 + cylHeight*cyl.scaleY + baseHeight*2,2* body1SizeZ/3);
	glScaled(0.75, 1, 1.0 / 3);
	body1.SetColor(5);
	if (bWireFrame)
		body1.DrawWireframe();
	else{
		setupMaterial(ambient, diffuse, specular, shiness);
		body1.Draw();
	}
		

	glPopMatrix();
}

void drawBody3(){
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(-0.25 * body1SizeX, body1SizeY * 2 + cylHeight*cyl.scaleY + baseHeight*2, -2 * body1SizeZ / 3);
	glScaled(0.75, 1, 1.0 / 3);

	body1.SetColor(5);
	if (bWireFrame)
		body1.DrawWireframe();
	else{
		setupMaterial(ambient, diffuse, specular, shiness);
		body1.Draw();
	}


	glPopMatrix();
}

void drawBody4(){
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	
	glTranslated(0.75 * body1SizeX, body1SizeY * 3 + cylHeight*cyl.scaleY + baseHeight*2,0);
	glScaled(0.25, 2.0, 1);
	body1.SetColor(6);
	if (bWireFrame)
		body1.DrawWireframe();
	else{
		setupMaterial(ambient, diffuse, specular, shiness);
		body1.Draw();
	}
	glPopMatrix();
}

void drawArm1(){
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 1, -1);
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0); 
	glTranslated(1.25 +1, 0.1 + body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight*2, 0);
	glRotatef(cottar1.rotateY, 0, 1, 0);
	glTranslated(-1, 0, 0);

	if (bWireFrame)
		arm1.DrawWireframe();
	else{
		
		setupMaterial(arm1ambient, diffuse, specular, shiness);

		arm1.Draw();
			
	}
	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawFrame1(){
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(1.25 + 1, 0.1 + body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight * 2, 0);
	glRotatef(cottar1.rotateY, 0, 1, 0);
	glTranslated(-1, 0, 0);
	if (selection == 1){
		GLfloat frameambient[] = { 1.0, 1.0, 1.0, 1.0f };
		glTranslated(0.15, 0, 0);
		setupMaterial(frameambient, diffuse, specular, shiness);
		frame1.DrawWireframe();
	}
	glPopMatrix();
}

void drawArm2(){
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS,2, -1);
	float alpha = asin(2.0*sin(cottar1.rotateY*PI / 180) / 3)*180.0 / PI;
	float beta = 180.0 - alpha - cottar1.rotateY;
	float x = sqrt(pow(2.0, 2) + pow(3.0, 2) - 2.0*2.0*3.0*cos(beta*PI / 180));
	glPushMatrix(); 
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(2.25 - x, 0, 0);
	glRotatef(-alpha,0,1,0);
	glTranslated(1.5, 0.3 + body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight * 2, 0);

	if (bWireFrame)
		arm2.DrawWireframe();
	else{
		setupMaterial(arm2ambient, diffuse, specular, shiness);
		arm2.Draw();
	}	
	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawFrame2(){
	float alpha = asin(2.0*sin(cottar1.rotateY*PI / 180) / 3)*180.0 / PI;
	float beta = 180.0 - alpha - cottar1.rotateY;
	float x = sqrt(pow(2.0, 2) + pow(3.0, 2) - 2.0*2.0*3.0*cos(beta*PI / 180));
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(2.25 - x, 0, 0);
	glRotatef(-alpha, 0, 1, 0);
	glTranslated(1.5, 0.3 + body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight * 2, 0);
	if (selection == 2){
		GLfloat frameambient[] = { 1.0, 1.0, 1.0, 1.0f };
		setupMaterial(frameambient, diffuse, specular, shiness);
		frame2.DrawWireframe();
	}
	glPopMatrix();
}

void drawCottar1(){
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 1, -1);
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(body1SizeX *0.375+1, body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight*2, 0);
	glRotatef(cottar1.rotateY, 0, 1, 0);

	if (bWireFrame)
		cottar1.DrawWireframe();
	else{
		GLfloat ambientcottar[4] = { 1,0,0, 0 };
		setupMaterial(ambientcottar, diffuse, specular, shiness);
		cottar1.Draw();
	}
		

	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawSlide(){
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 3, -1);
	float transy = 0.1 + body1SizeY * 3 + cylHeight*cyl.scaleY + baseHeight * 2;
	float alpha = asin(2.0*sin(cottar1.rotateY*PI/180)/3)*180.0/PI;
	float beta = 180.0 - alpha - cottar1.rotateY;
	float x = sqrt(pow(2.0, 2) + pow(3.0, 2) - 2.0*2.0*3.0*cos(beta*PI / 180));
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(2.25 - x , transy, 0);

	body1.SetColor(5);
	if (bWireFrame)
		slide.DrawWireframe();
	else{
		setupMaterial(slideambient, diffuse, specular, shiness);
		slide.Draw();
	}
	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawFrame3(){
	float transy = 0.1 + body1SizeY * 3 + cylHeight*cyl.scaleY + baseHeight * 2;
	float alpha = asin(2.0*sin(cottar1.rotateY*PI / 180) / 3)*180.0 / PI;
	float beta = 180.0 - alpha - cottar1.rotateY;
	float x = sqrt(pow(2.0, 2) + pow(3.0, 2) - 2.0*2.0*3.0*cos(beta*PI / 180));
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	glTranslated(2.25 - x, transy, 0);
	if (selection == 3){
		GLfloat frameambient[] = { 1.0, 1.0, 1.0, 1.0f };
		setupMaterial(frameambient, diffuse, specular, shiness);
		frame3.DrawWireframe();
	}
	glPopMatrix();
}

void drawCottar2(){
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 2, -1);
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);
	
	glTranslated(1.25 + 1, 0.2+ body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight * 2, 0);
	glRotatef(cottar1.rotateY, 0, 1, 0);
	glTranslated(-2, 0, 0);

	if (bWireFrame)
		cottar2.DrawWireframe();
	else{
		GLfloat ambientcottar[4] = { 1, 0, 0, 0 };
		setupMaterial(ambientcottar, diffuse, specular, shiness);
		cottar2.Draw();
	}
		

	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawCottar3(){
	glEnable(GL_STENCIL_TEST);
	glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
	glStencilFunc(GL_ALWAYS, 2, -1);
	float transy = 0.2 + body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight * 2;
	float alpha = asin(2.0*sin(cottar1.rotateY*PI / 180) / 3)*180.0 / PI;
	float beta = 180.0 - alpha - cottar1.rotateY;
	float x = sqrt(pow(2.0, 2) + pow(3.0, 2) - 2.0*2.0*3.0*cos(beta*PI / 180));
	glPushMatrix();
	glRotatef(base.rotateY, 0, 1, 0);

	//glTranslated(body1SizeX *0.375 + 1, 0.2 + body1SizeY * 5 + cylHeight*cyl.scaleY + baseHeight * 2, 0);
	glTranslated(2.25 - x, transy, 0);

	if (bWireFrame)
		cottar2.DrawWireframe();
	else{
		GLfloat ambientcottar[4] = { 1, 0, 0, 0 };
		setupMaterial(ambientcottar, diffuse, specular, shiness);
		cottar2.Draw();
	}
		

	glPopMatrix();
	glDisable(GL_STENCIL_TEST);
}

void drawCyl(){
	glPushMatrix();
	glTranslated(0, cyl.slideY + baseHeight, 0);
	glScalef(cyl.scaleX, cyl.scaleY, cyl.scaleZ);
	glRotatef(base.rotateY, 0, 1, 0);
	if (bWireFrame)
		cyl.DrawWireframe();
	else{
		GLfloat cylambient[] = { 0.0, 0.0, 1.0, 1.0f };
		setupMaterial(cylambient, diffuse, specular, shiness);
		cyl.Draw();
	}
		
	glPopMatrix();
}

void drawFloor(){
	//glTranslated(0, 0.001, 0);
	loadTextures();
	float fXWidth = 50;
	float fZWidth = 50;
	glEnable(GL_TEXTURE_2D);
	
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glBindTexture(GL_TEXTURE_2D, floorTex.texID);
	glColor4f(1, 1, 1, 0.5);
	
	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 0.0);
	glVertex3f(-fXWidth / 2, -0.001, fZWidth / 2);

	glTexCoord2f(1.0, 0);
	glVertex3f(fXWidth / 2, -0.001, fZWidth / 2);

	glTexCoord2f(1.0, 1.0);
	glVertex3f(fXWidth / 2, -0.001, -fZWidth / 2);

	glTexCoord2f(0.0, 1.0);
	glVertex3f(-fXWidth / 2, -0.001, -fZWidth / 2);
	glEnd();
	
	glDisable(GL_TEXTURE_2D);
	
	glPopMatrix();
	glEnable(GL_LIGHTING);	
	
}

void drawShadow(bool on_off){
	glDisable(GL_LIGHTING);
	
	glPushMatrix();

	GLfloat round[4] = { 0, 1, 0, 0 };
	GLfloat light[4] = { 6, 6, 6, 0 };
	glColor4f(0.5, 0.5, 0.5, 1);
	myShadowMatrix(round, light);

	drawBase();
	drawCyl();
	drawBody1();
	drawBody2();
	drawBody3();
	drawBody4();
	drawCottar1();
	drawCottar2();
	drawCottar3();
	drawArm1();
	drawArm2();
	drawSlide();
	glPopMatrix();
	glEnable(GL_LIGHTING);
	if (on_off){
		glDisable(GL_LIGHTING);

		glPushMatrix();

		GLfloat round1[4] = { 0, 1, 0, 0 };
		GLfloat light1[4] = { 6.0f, 6.0f, -6.0f, 0.0f };
		glColor4f(0.5, 0.5, 0.5, 1);
		myShadowMatrix(round1, light1);

		drawBase();
		drawCyl();
		drawBody1();
		drawBody2();
		drawBody3();
		drawBody4();
		drawCottar1();
		drawCottar2();
		drawCottar3();
		drawArm1();
		drawArm2();
		drawSlide();
		glPopMatrix();
		glEnable(GL_LIGHTING);
	}
}

void light1(bool on_off){
	GLfloat	lightDiffuse1[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat	lightSpecular1[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat	lightAmbient1[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat light_position2[] = { 6.0f, 6.0f, -6.0f, 0.0f };
	glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular1);
	if (on_off){
		glEnable(GL_LIGHT1);
	}
	else {
		glDisable(GL_LIGHT1);

	}
	
}

void light0(){
	GLfloat	lightDiffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat	lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat	lightAmbient[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat light_position1[] = { 6.0f, 6.0f, 6.0f, 0.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}

void alldraw(){
	light0();
	light1(on_off);
	
	drawBase();
	drawCyl();
	drawBody1();
	drawBody2();
	drawBody3();
	drawBody4();
	drawCottar1();
	drawCottar2();
	drawCottar3();
	
	drawArm1();
	drawFrame1();
	drawArm2();
	drawFrame2();
	drawSlide();
	drawFrame3();

}

void myDisplay(){
	
	glMatrixMode(GL_MODELVIEW);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glViewport(0, 0, screenWidth, screenHeight);

	drawShadow(on_off);
	setupMaterial(ambient, diffuse, specular, shiness);
	alldraw();
	
	glPushMatrix();
	glScalef(1.0, -1.0, 1.0);
	alldraw();
	glPopMatrix();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	drawFloor();
	glDisable(GL_BLEND);
	
	glutSwapBuffers();
	glFlush();
}

void myInit(){
	float	fHalfSize = 4;
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glColor3f(1.0f, 1.0f, 1.0f);

	glFrontFace(GL_CCW);
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glEnable(GL_NORMALIZE);
	glLoadIdentity();
	glOrtho(-fHalfSize, fHalfSize, -fHalfSize, fHalfSize, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);

}

int _tmain(int argc, _TCHAR* argv[]){

	cout << "1, 2: Rotate the base" << endl;
	cout << "3, 4: Rotate the arm" << endl;
	cout << "L, l: Cylinder up" << endl;
	cout << "X, x: Cylinder down" << endl;
	cout << "W, w: Switch between wireframe and solid mode" << endl;

	glutInit(&argc, (char**)argv); //initialize the tool kit
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);//set the display mode
	glutInitWindowSize(screenWidth, screenHeight); //set window size
	glutInitWindowPosition(100, 100); // set window position on screen
	glutCreateWindow("ASM1"); // open the screen window
	
	base.CreateCylinder(20, baseHeight, baseRadius);
	base.CalculateFacesNorm();
	base.SetColor(7);
	base.slideY = baseHeight / 2.0;
	
	
	
	cyl.CreateCylinder(20, cylHeight, cylRadius);
	cyl.SetColor(9);
	cyl.scaleX = 1.0;
	cyl.scaleY = 1.0;
	cyl.scaleZ = 1.0;
	cyl.CalculateFacesNorm();


	body1.CreateCuboid(body1SizeX, body1SizeY, body1SizeZ);
	body1.SetColor(2);
	body1.CalculateFacesNorm();


	arm1.CreateSlipper(0.6, 0.3, 2.0, 0.2);
	arm1.SetColor(10);
	arm1.CalculateFacesNorm();
	frame1.CreateCuboid(1.45, 0.1, 0.6);
	
	arm2.CreateSlipper(0.3, 0.3, 3, 0.2);
	arm2.SetColor(8);
	arm2.CalculateFacesNorm();
	frame2.CreateCuboid(1.8, 0.1, 0.3);

	cottar1.CreateCylinder(20, 0.21, cylRadius);
	cottar1.SetColor(11);
	cottar1.CalculateFacesNorm();


	cottar2.CreateCylinder(20, 0.21, cylRadius / 2);
	cottar2.SetColor(11);
	cottar2.CalculateFacesNorm();


	cottar3.CreateCylinder(20, 0.21, cylRadius / 2);
	cottar2.SetColor(11);
	cottar3.CalculateFacesNorm();


	slide.CreateCuboid(body1SizeX/12, body1SizeY * 2+0.1, body1SizeZ/3);
	slide.SetColor(12);
	slide.CalculateFacesNorm();
	frame3.CreateCuboid(body1SizeX / 12, body1SizeY * 2 + 0.1, body1SizeZ / 3);

	myInit();
	glutSpecialFunc(keyboardFunc);
	glutKeyboardFunc(myKeyboard);
	glutMouseFunc(onMouse);
	glutDisplayFunc(myDisplay);

	cam.set(Point3(-8, 4, 8), Point3(0, 1, 0), Vector3(0, 1, 0));
	cam.setShape(30.0f, 64.0f / 48.0f, 0.5f, 500.0f);
	
	glutMainLoop();
	return 0;
}