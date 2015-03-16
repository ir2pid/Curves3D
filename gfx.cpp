
//Include OpenGL header files, so that we can use OpenGL
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <string.h>
#include <GL/glui.h>
#include <iostream>
#include <stdlib.h> //Needed for "exit" function
#include <math.h>   //Needed for "pow" function


#pragma comment( linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"" )  //to hide inactive DOS console wondow


int x=0;
using namespace std;

#define SELECT 1
#define DUPLICATE 2
#define DELETE 3
#define DELETE_ALL 4
#define REFRESH 5


int PARAMETRIC_BEIZER =0;
int SUBDIVISION_BEIZER =1;

int promptHUDPos = -50;
int promptHUDShow = 0;
int autoUpdate = 1;
int sketchEffect=0;
int showKey = 1;
int showKeyPos = -50;
int showGuide = 1;
int introID=1;
int showIntro =1;
int display_3D = 1;
int showPoints = 0;
int c1Continuity = 1;
int radiogroup_item_id=0;
int RADIO=0;
int segments=1;
int strength = 100;
double tStrength=0.0;
int x_pixel = 0;
int y_pixel = 0;
int snap = 0;
int x_cord = 0;
int y_cord = 0;
bool fullscreen = 0;
bool select_on = 0;
bool duplicate_on = 0;
bool toggle = 0;
int point_id = 0;

int   main_window;
float promptHUDColor=0.0f;
float _angle = 30.0f;
float _cameraAngle = 0.0f;
GLUI *glui;
GLUI_Panel       *hudOptions,*curveOptions, *cordDisplay;
GLUI_RadioGroup  *radio_group;
GLUI_Spinner *spinner,*spinner2,*spinnerIntroID;

//---------------------------------------
class point
{
public:

	double x,y;
	point()
	{
		x=0;
		y=0;
	}
	point(double a,double b)
	{
		x=a;
		y=b;
	}
	void set(point m)
	{
		x=m.x;
		y=m.y;
	}

};

class cubeMatrix
{
private:
	double Mx[4][4];
	double Gx[4][1]; //holds px
	double My[4][4];
	double Gy[4][1]; //holds py

	double B[1][4]; //holds t [0,1]
	double tempx[4][1];
	double tempy[4][1];
	double tempx2[1][4];
	double tempy2[1][4];
	int sx,sy;
	point out;
	
public:
	cubeMatrix()
	{
		Mx[0][0]=-1;	Mx[0][1]=3;		Mx[0][2]=-3;	Mx[0][3]=1;
		Mx[1][0]=3;		Mx[1][1]=-6;	Mx[1][2]=3;		Mx[1][3]=0;
		Mx[2][0]=-3;	Mx[2][1]=0;		Mx[2][2]=3;		Mx[2][3]=0;
		Mx[3][0]=1;		Mx[3][1]=4;		Mx[3][2]=1;		Mx[3][3]=0;

		My[0][0]=-1;	My[0][1]=3;		My[0][2]=-3;	My[0][3]=1;
		My[1][0]=3;		My[1][1]=-6;	My[1][2]=3;		My[1][3]=0;
		My[2][0]=-3;	My[2][1]=0;		My[2][2]=3;		My[2][3]=0;
		My[3][0]=1;		My[3][1]=4;		My[3][2]=1;		My[3][3]=0;
	}
	
	void initG(vector<point> temp)
	{
		for(int i=0;i<temp.size();i++)
		{
			Gx[i][0]=temp[i].x;
			Gy[i][0]=temp[i].y;
		}
	}
	
	void initB(double temp)
	{
		for(int i=2;i>=0;i--)
		{
			B[0][i]=temp;
			temp=temp*B[0][2];
		}
			B[0][3]=1.0f;
			
			
	}

	point matrixMul()
	{

		tempx[0][0]=(Mx[0][0]*Gx[0][0])+(Mx[0][1]*Gx[1][0])+(Mx[0][2]*Gx[2][0])+(Mx[0][3]*Gx[3][0]);
		tempx[1][0]=(Mx[1][0]*Gx[0][0])+(Mx[1][1]*Gx[1][0])+(Mx[1][2]*Gx[2][0])+(Mx[1][3]*Gx[3][0]);
		tempx[2][0]=(Mx[2][0]*Gx[0][0])+(Mx[2][1]*Gx[1][0])+(Mx[2][2]*Gx[2][0])+(Mx[2][3]*Gx[3][0]);
		tempx[3][0]=(Mx[3][0]*Gx[0][0])+(Mx[3][1]*Gx[1][0])+(Mx[3][2]*Gx[2][0])+(Mx[3][3]*Gx[3][0]);

		tempy[0][0]=(My[0][0]*Gy[0][0])+(My[0][1]*Gy[1][0])+(My[0][2]*Gy[2][0])+(My[0][3]*Gy[3][0]);
		tempy[1][0]=(My[1][0]*Gy[0][0])+(My[1][1]*Gy[1][0])+(My[1][2]*Gy[2][0])+(My[1][3]*Gy[3][0]);
		tempy[2][0]=(My[2][0]*Gy[0][0])+(My[2][1]*Gy[1][0])+(My[2][2]*Gy[2][0])+(My[2][3]*Gy[3][0]);
		tempy[3][0]=(My[3][0]*Gy[0][0])+(My[3][1]*Gy[1][0])+(My[3][2]*Gy[2][0])+(My[3][3]*Gy[3][0]);

		out.x = ( (B[0][0]*tempx[0][0])+(B[0][1]*tempx[1][0])+(B[0][2]*tempx[2][0])+(B[0][3]*tempx[3][0]) )/6.0;
		out.y = ( (B[0][0]*tempy[0][0])+(B[0][1]*tempy[1][0])+(B[0][2]*tempy[2][0])+(B[0][3]*tempy[3][0]) )/6.0;

		return out;

	}

};

class quadMatrix
{
private:
	double Mx[3][3];
	double Gx[3][1]; //holds px
	double My[3][3];
	double Gy[3][1]; //holds py

	double B[1][3]; //holds t [0,1]
	double tempx[3][1];
	double tempy[3][1];
	double tempx2[1][3];
	double tempy2[1][3];
	int sx,sy;
	point out;
	
public:
	quadMatrix()
	{
		Mx[0][0]=1;		Mx[0][1]=-2;	Mx[0][2]=1;	
		Mx[1][0]=-2;	Mx[1][1]=2;		Mx[1][2]=0;		
		Mx[2][0]=1;		Mx[2][1]=1;		Mx[2][2]=0;		
		

		My[0][0]=1;		My[0][1]=-2;	My[0][2]=1;	
		My[1][0]=-2;	My[1][1]=2;		My[1][2]=0;		
		My[2][0]=1;		My[2][1]=1;		My[2][2]=0;	
	}
	
	void initG(vector<point> temp)
	{
		for(int i=0;i<temp.size();i++)
		{
			Gx[i][0]=temp[i].x;
			Gy[i][0]=temp[i].y;
		}
	}
	
	void initB(double temp)
	{
		for(int i=1;i>=0;i--)
		{
			B[0][i]=temp;
			temp=temp*B[0][1];
		}
			B[0][2]=1.0f;
			
			
	}

	point matrixMul()
	{

		tempx[0][0]=(Mx[0][0]*Gx[0][0])+(Mx[0][1]*Gx[1][0])+(Mx[0][2]*Gx[2][0]);
		tempx[1][0]=(Mx[1][0]*Gx[0][0])+(Mx[1][1]*Gx[1][0])+(Mx[1][2]*Gx[2][0]);
		tempx[2][0]=(Mx[2][0]*Gx[0][0])+(Mx[2][1]*Gx[1][0])+(Mx[2][2]*Gx[2][0]);
		

		tempy[0][0]=(My[0][0]*Gy[0][0])+(My[0][1]*Gy[1][0])+(My[0][2]*Gy[2][0]);
		tempy[1][0]=(My[1][0]*Gy[0][0])+(My[1][1]*Gy[1][0])+(My[1][2]*Gy[2][0]);
		tempy[2][0]=(My[2][0]*Gy[0][0])+(My[2][1]*Gy[1][0])+(My[2][2]*Gy[2][0]);
		

		out.x = ( (B[0][0]*tempx[0][0])+(B[0][1]*tempx[1][0])+(B[0][2]*tempx[2][0]) )/2.0;
		out.y = ( (B[0][0]*tempy[0][0])+(B[0][1]*tempy[1][0])+(B[0][2]*tempy[2][0]) )/2.0;

		return out;

	}

};



//---------------------------------------

vector<int> ivx,ivy;
vector<point> curve;


void orthographic_projection()
{
	// switch to projection mode
	glMatrixMode(GL_PROJECTION);
	// save previous matrix which contains the 
	//settings for the perspective projection
	glPushMatrix();
	// reset matrix
	glLoadIdentity();
	// set a 2D orthographic projection
	gluOrtho2D(0, x_pixel , 0, y_pixel);					// w, h set by window_parameters
	// invert the y axis, down is positive
	glScalef(1, -1, 1);
	// mover the origin from the bottom left corner
	// to the upper left corner
	glTranslatef(0, -y_pixel, 0);
	glMatrixMode(GL_MODELVIEW);

}

void perspective_projection()
{
	glMatrixMode(GL_PROJECTION);

	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);

}



void renderBitmapCharacter(float x, float y, float z, void *font,char *string)
{
  
  char *c;
  glRasterPos2f(x, y);
  for (c=string; *c != '\0'; c++) {
    glutBitmapCharacter(font, *c);
  }
}

void renderCircle(int x0, int y0, int radius)
{
		int x1=0;
		int y1=radius;
		int d=3-(2*radius);
		glBegin(GL_POINTS);
		while(x1<=y1)
			{
			glVertex2f(x0+x1,y0+y1);
			glVertex2f(x0-x1,y0+y1);
			glVertex2f(x0+x1,y0-y1);
			glVertex2f(x0-x1,y0-y1);
			glVertex2f(x0+y1,y0+x1);
			glVertex2f(x0-y1,y0+x1);
			glVertex2f(x0+y1,y0-x1);
			glVertex2f(x0-y1,y0-x1);
			if(d<0)
			d+=4*x1+6;
			else
				{
				d+=4*(x1-y1)+10;
				y1--;
				}
			x1++;
			}
		glEnd();
}


void promptHUD(char *c)
{

	char tmp[100];

	sprintf(tmp,"%s",c);
	if( toggle ==0 )
	{
		
		if(promptHUDColor > 1.0f)
		{
			toggle = 1;
		}
		else
		{
			promptHUDColor+=0.05;
		}
	}
	if( toggle==1 )
	{

		if(promptHUDColor < 0.0f)
		{
			toggle=2;
		}
		else
		{
			promptHUDColor-=0.05;
		}
	}

			glColor3f(0.0f,promptHUDColor,0.0f);
			renderBitmapCharacter(420,17,0,GLUT_BITMAP_8_BY_13,tmp);
			glBegin(GL_POINTS);
			glVertex2f(415.0f,15.0f);
			glEnd();
			renderCircle(415,15,2);


}

/*
void promptHUD(char *c)
{

	char tmp[100];
	int x;
	sprintf(tmp,"%s",c);

	if(promptHUDPos<9)
	{
		x=1;
	}
	if(promptHUDPos>8 )
	{
		x=-50;
	}
		promptHUDPos+=x;

			glColor3f(0.0f,0.5f,1.0f);
			renderBitmapCharacter(400,17+promptHUDPos,0,GLUT_BITMAP_8_BY_13,tmp);
			glBegin(GL_POINTS);
			glVertex2f(400.0f,15.0f+promptHUDPos);
			glEnd();
			renderCircle(400,15+promptHUDPos,2);


}*/

// generate cubic bspline
void genCubeBSplineCurve()
{
	if(ivx.size()>3)
	{
		curve.clear();

			for(int i=0;i<ivx.size()-1;i++)
			{
				if((i+3) <= ivx.size()-1)
				{
					for(double t = 0.0; t < 1.0; t += tStrength) 
					{
						point p;
						vector<point> vp;
						cubeMatrix mat;
						mat.initB(t);
						for(int j=i;j<=i+3;j++)
						{
							p.x=ivx[j];
							p.y=ivy[j];
							vp.push_back(p);
						}
						mat.initG(vp);
						p=mat.matrixMul();
						curve.push_back(p);
					}
				}
			}
	}
}


// generate quad bspline
void genQuadBSplineCurve()
{
	if(ivx.size()>2)
	{
		curve.clear();

			for(int i=0;i<ivx.size()-1;i++)
			{
				if((i+2) <= ivx.size()-1)
				{
					for(double t = 0.0; t < 1.0; t += tStrength) 
					{
						point p;
						vector<point> vp;
						quadMatrix mat;
						mat.initB(t);
						for(int j=i;j<=i+2;j++)
						{
							p.x=ivx[j];
							p.y=ivy[j];
							vp.push_back(p);
						}
						mat.initG(vp);
						p=mat.matrixMul();
						curve.push_back(p);
					}
				}
			}
	}
}


int c(double n,double r)
{
	//c(n,r) combination function
	if(r == 0 || n==r)
	{
		return 1.0;
	}
	else
	{	
		double x = 1.0;
		for(double i = r; i > 0; i--)
		{
			x = x * n;
			n--;
		}
		for(double i = r; i > 0; i--)
		{
			x = x / i;
		}
		return x;
	}
}


//generate the parametric curve
void genParametricCurve()
{
	if(ivx.size() > 1)
		{
			curve.clear();

			for(double t = 0.0;t < 1.0; t += tStrength) 
			{
				int n = ivx.size()-1;
				int r = 0;		
				point p;
				p.x = 0;
				p.y = 0;

				for(int i=0;i<ivx.size();i++)
				{
					double x,y,z;
					x = c((double)n,(double)r);
					y = pow(t,r);
					z = pow((1-t),n-r);

						p.x = p.x + x * y * z * (double)ivx[i] ; 
						p.y = p.y + x * y * z * (double)ivy[i] ;


					r++;
				}

			curve.push_back(p);
		}
	}
}



void genCubeCurve(vector<point> input)
{

	
	for(double t = 0.0;t < 1.0; t += tStrength) 
	{
		point p;
		p.x = ( (pow((1 - t), 3) * (float)input[0].x) + (3 * t * (pow((1 - t), 2)) * (float)input[1].x) + (3 * (1-t) * pow(t, 2) * (float)input[2].x) + pow(t,3) * input[3].x);
		p.y = ( (pow((1 - t), 3) * (float)input[0].y) + (3 * t * (pow((1 - t), 2)) * (float)input[1].y) + (3 * (1-t) * pow(t, 2) * (float)input[2].y) + pow(t,3) * input[3].y);
		curve.push_back(p);
    }
}

void genC1CubeCurve()
{
	curve.clear();
	if(ivx.size()>3)
	{
		vector<point> input;
		point p;
		for(int i=0;i+3<ivx.size();i=i+3)
		{
			for(int j=i;j<i+4;j++)
			{
				p.x=ivx[j];
				p.y=ivy[j];

				//enforcing c1 continuity
	
				if(RADIO==5 && j>3 && (j-1)%3==0)
				{	//q1 = pn - (pn-1) + q0
					p.x=(2*ivx[j-1]) - ivx[j-2];
					p.y=(2*ivy[j-1]) - ivy[j-2];
					ivx[j]=p.x;
					ivy[j]=p.y;
				}
				input.push_back(p);
			}
			genCubeCurve(input);
			input.clear();
		}
	}
}

void genQuadCurve(vector<point> input)
{
	

	for(double t = 0.0;t < 1.0; t += tStrength) 
	{
		point p;
		p.x = pow((1 - t), 2) * (float)input[0].x + 2 * t * (1 -t) * (float)input[1].x + pow(t, 2) * (float)input[2].x;
		p.y = pow((1 - t), 2) * (float)input[0].y + 2 * t * (1 -t) * (float)input[1].y + pow(t, 2) * (float)input[2].y;
		curve.push_back(p);
    }
}

void genC1QuadCurve()
{
	curve.clear();
	if(ivx.size()>2)
	{
		vector<point> input;
		point p;
		for(int i=0;i+2<ivx.size();i=i+2)
		{
			for(int j=i;j<i+3;j++)
			{
				p.x=ivx[j];
				p.y=ivy[j];

				//enforcing c1 continuity
		
				if(RADIO==7 && j>2 && (j-1)%2==0)
				{	//q1 = pn - (pn-1) + q0
					p.x=(2*ivx[j-1]) - ivx[j-2];
					p.y=(2*ivy[j-1]) - ivy[j-2];
					ivx[j]=p.x;
					ivy[j]=p.y;
				}
				input.push_back(p);
			}
			genQuadCurve(input);
			input.clear();
		}
	}
}


void OneSubdivide(vector<point> p ,vector<point> poly1,vector<point> poly2,double u )
{
	
	curve.clear();
	if(p.size()!=0)
	{	
		vector<point> output;
		if (p.size()==1)
		{		
			//
		
			point p0;
			p0.x=p[0].x;
			p0.y=p[0].y;
			//store in output vector poly1.p0.poly2
			output.insert(output.end(),poly1.begin(),poly1.end());
			output.push_back(p0);
			output.insert(output.end(),poly2.begin(),poly2.end());
			curve.reserve(output.size());
			curve.insert(curve.begin(),output.begin(),output.end());
			//return output;
		}
		else
		{
			
			point p0,pn;

			//poly1=poly1.p0
			p0.x=p[0].x;
			p0.y=p[0].y;
			poly1.push_back(p0);

			//poly2=pn.poly2
			pn.x=p[(p.size()-1)].x;
			pn.y=p[(p.size()-1)].y;
			poly2.insert(poly2.begin(),pn);


			int size=p.size();
			double varX,varY;
			for (int i=0;i <=(size-2);i++)
			{
				point temp;
				temp.x=p[i].x+ u*(p[i+1].x-p[i].x);
				temp.y=p[i].y+ u*(p[i+1].y-p[i].y);
				output.push_back(temp);
			}
		
			//call OneSubdivide again
			OneSubdivide(output,poly1,poly2,u);

		}


	}
}
//genereate curve by subdivisionmethod
void genSubDivisionCurve(vector<point> in ,int m,double u)//this is basically One subdivide
{
	if(ivx.size() > 1)
	{
		//tryingto insert at the beginnig of the vector
		vector<int> a,b,p3x,p3y,p4x,p4y;
		vector<point> p1,p2,input,output;


		input.insert(input.begin(),in.begin(),in.end());

		if (m==1)
		{
				OneSubdivide(input,p1,p2,0.5);
		}
		else
		{
			vector<point> input1,input2;
			OneSubdivide(input,p1,p2,0.5);
			int siz=curve.size();
			int n=(siz/2) + 1;

			int i;
			for (i=0;i<n;i++)
			{
				point temp(curve[i].x,curve[i].y);
				input1.push_back(temp);
			}

	
			for (--i;i<siz;i++)
			{
				point temp(curve[i].x,curve[i].y);
				input2.push_back(temp);
			}
			vector<point> final;
			final.clear();
			genSubDivisionCurve(input2,m-1,0.5);
			final.insert(final.begin(),curve.begin(),curve.end());
			genSubDivisionCurve(input1,m-1,0.5);
			final.insert(final.begin(),curve.begin(),curve.end());
			curve.clear();
			curve.insert(curve.begin(),final.begin(),final.end());			
		}
	}

}

void genBSplinesSubdivide()
{
	point p,p1,p2,p3;
	int size;
	
	for (int m=0;m<ivx.size();m++)
	{
		p.x=ivx[m];p.y=ivy[m];
		curve.push_back(p);
	}

	for(int i=0;i<segments;i++)
	{
		size=curve.size();
		//std::cout<<" before "<<size<<endl;
		if(size>=3)
		{
			for(int j=0;j<=size-3;j++)
			{
				if (j==0)
				{
				p1.x=(curve[j].x + curve[j+1].x)/2;
				p1.y=(curve[j].y + curve[j+1].y)/2;curve.push_back(p1);}

				p2.x=(curve[j].x/8)+(3*curve[j+1].x/4)+(curve[j+2].x/8);
				p2.y=(curve[j].y/8)+(3*curve[j+1].y/4)+(curve[j+2].y/8);curve.push_back(p2);
				p3.x=(curve[j+1].x + curve[j+2].x)/2;
				p3.y=(curve[j+1].y + curve[j+2].y)/2;curve.push_back(p3);


			}

			curve.erase(curve.begin(),curve.begin()+size);
		}
	 }
}

//generate genereic curve
void genCurve()
{	
	/// all this happens when we do parametric
	switch(RADIO)
	{
		case 0 :
			{
				curve.clear();
				genQuadBSplineCurve();
				break;				
			}		
		case 1 :
			{
				curve.clear();
				genCubeBSplineCurve();
				break;				
			}	
		case 2:
			{
				curve.clear();
				genBSplinesSubdivide();
				break;
			}
		case 3 :
			{	
				curve.clear();
				genParametricCurve();
				break ;
			}
		case 4 :
			{	//c0 checked in genC1CubeCurve()
				curve.clear();
				genC1CubeCurve();
				break ;
			}

		case 5 :
			{	//c1 checked in genC1CubeCurve()
				curve.clear();
				genC1CubeCurve();
				break ;
			}


		case 6 :
			{	//c1 checked in genC1QuadCurve()
				curve.clear();
				genC1QuadCurve();
				break ;
			}
		case 7 :
			{	//c1 checked in genC1QuadCurve()
				curve.clear();
				genC1QuadCurve();
				break ;
			}
		case 8 :
			{
				curve.clear();
				vector<point> input;
				for(int i=0;i<ivx.size();i++)	// to store ivx,ivy in vector<point> input
					{	
						//constructor ovrldin
						point temp(ivx[i],ivy[i]);
						input.push_back(temp);
					}
				genSubDivisionCurve(input,segments,0.5);
				break;
			}

		default:break;

	}
	
}


void drawScene() {
	

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(display_3D)
	{
		
		glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
		glLoadIdentity(); //Reset the drawing perspective
		glRotatef(-_cameraAngle, 0.0f, 0.0f, 0.0f); //Rotate the camera
		glTranslatef(0.0f, 0.0f, -5.0f); //Move forward 5 units
		
		
			float rndA=((rand()% 100))/10 +0.01;
			float rndB=((rand()% 100))/10 +0.01;
			float rndC=((rand()% 100))/10 +0.01;

		glPushMatrix(); //Save the current state of transformations
			glTranslatef(-1.0f, 0.0f, 0.0f); //Move to the center of the pentagon
			glRotatef(_angle, rndC/50, 1.0f, rndA/50); //Rotate about the y-axis

				glColor3f(rndA/10,rndB/10+0.7f,rndC/10+ (_angle/360) );
			//glColor3f(0.0f,1.0f,0.0f);
			if(rndC/10 <0.5)
			glScalef(rndA,rndA,rndA);
			if(introID==1)
			{				
			//glutWireDodecahedron();
			glutWireIcosahedron();
			glPushMatrix();
			glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
			glutWireIcosahedron();
			glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
			glutWireIcosahedron();
			glPopMatrix();
			//glutWireOctahedron();
			glutWireTorus(0.1f, 2.0f , 5,50);
			glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
			glutWireTorus(0.1f, 2.0f , 5,50);
			glutWireTetrahedron();
			}
			if(introID==2)
			{				
			glutWireTorus(0.1f, 2.0f , 5,50);
			glutWireTeapot(1.0f);
			}
			if(introID==3)
			{
				
				glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
				glutWireSphere(1.0f,25,25);
				glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
				glutWireTorus(0.1f, 2.0f , 5,50);
				glRotatef(90.0f, 0.0f, 1.0f, 0.0f);
				glutWireTorus(0.1f, 2.0f , 5,50);
				glRotatef(90.0f, 1.0f, 0.0f, 0.0f);
				glutWireTorus(0.1f, 2.0f , 5,50);
			}
		

		
		glPopMatrix();

		glPushMatrix();
				glColor3f(rndA/10,rndB/10+0.7f,rndC/10+ (_angle/360) );
		if(rndC/10 <0.8)
		renderBitmapCharacter(-1.16-(rndA/99),1.7-(rndB/99),0,GLUT_BITMAP_9_BY_15,"Lab 1");
		else
		renderBitmapCharacter(-1.16,1.7,0,GLUT_BITMAP_9_BY_15,"Lab 1");
		//glColor3f(0.4f,0.0f,0.0f);
		//renderBitmapCharacter(-0.12,1.7,0,GLUT_BITMAP_8_BY_13,"Lab 1");
		glPopMatrix();
		glColor3f(1.0f,1.0f,1.0f);
		
			//motion blur
			float q = .90;
			glAccum(GL_MULT, q);
			glAccum(GL_ACCUM, 1-q);
			glAccum(GL_RETURN, 1.0);
			//motion blur
	}
	
	//setting 2d hud display and printing click locations
	orthographic_projection();
	glPushMatrix();
	glLoadIdentity();

	

	//draw points
	for(int i=0;i<(int)ivx.size();i++)
	{
		if( (abs( ivx[i]-x_cord) < 4)  && (abs( ivy[i]-y_cord) < 4 ) )
		{
			//glColor3f(0.0f,1.0f,0.0f);
			//renderCircle(ivx[i] ,ivy[i],4);
			glColor3f(1.0f,0.0f,0.0f);
		}
		else
		{
			glColor3f(1.0f,1.0f,1.0f);
		}

		glBegin(GL_POINTS);
			glVertex2f( (float)ivx[i] ,(float)ivy[i] );
		glEnd();


		if(snap && i == point_id)
		{
			renderCircle(ivx[i],ivy[i],5);
		}
		else
		{
			renderCircle(ivx[i],ivy[i],2);
		}

		glColor3f(1.0f,1.0f,1.0f);

	}



	if(curve.size()>0 )
	{

   //curve type based on RADIO
		switch(RADIO)
		{

		case 0:
			{	//quad b-spline
				glColor3f(0.0f,0.5f,1.0f);
				glBegin(GL_POINTS);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				}
				glEnd();
				glColor3f(1.0f,1.0f,1.0f);
				break;
			}	
		case 1:
			{	//cube b-spline
				glColor3f(0.0f,0.5f,1.0f);
				glBegin(GL_POINTS);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				}
				glEnd();
				glColor3f(1.0f,1.0f,1.0f);
				break;
			}	

			case 2:
			{	//bspline subdivision
				glColor3f(1.0f,0.3f,0.3f);
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);

				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				}
				glEnd();

				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);
				break;
			}

		case 3:
			{	
				//parametric
				glColor3f(1.0f,0.3f,1.0f);
				
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				
				}
				glEnd();
				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);
				break;
			}	
		case 4:
			{	
				//Cubic /w C1
				glColor3f(0.6f,0.6f,1.0f);
				
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				
				}
				glEnd();
				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);

				break;
			}	
		case 5:
			{	
				//Cubic /w C1
				glColor3f(0.6f,0.6f,1.0f);
				
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				
				}
				glEnd();
				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);

				break;
			}	
		
		case 6:
			{	
				//Quad /w C0
				glColor3f(1.0f,1.0f,0.0f);
				
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				
				}
				glEnd();
				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);

				break;
			}	

		case 7:
			{	
				//Quad /w C1
				glColor3f(1.0f,1.0f,0.0f);
				
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);
				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				
				}
				glEnd();
				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);

				break;
			}	


		case 8:
			{	//subdivision
				glColor3f(1.0f,0.5f,0.0f);
				glEnable(GL_LINE_STIPPLE);
				glLineStipple(4,0x5555);
				glBegin(GL_LINE_STRIP);

				for(int i=0;i<curve.size();i++)
				{
					glVertex2f(curve[i].x,curve[i].y);
				}
				glEnd();

				glDisable(GL_LINE_STIPPLE);
				glColor3f(1.0f,1.0f,1.0f);
				break;
			}
		default:break;



		}

	}

	// show control points
	if(showPoints)
	{
		glColor3f(0.0f,0.5f,1.0f);
		glBegin(GL_POINTS);
		for(int i=0;i<curve.size();i++)
		{
			glVertex2f(curve[i].x,curve[i].y);
			renderCircle(curve[i].x,curve[i].y,2);
		}
		glEnd();
		glColor3f(1.0f,1.0f,1.0f);
	}

	// show control guide
	if(showGuide)
	{

	//promptHUD("Hello");

		//guide lines
			glColor3f(0.0f,1.0f,0.0f);
			glEnable(GL_LINE_STIPPLE);
			glLineStipple(4,0x5555);
			glBegin(GL_LINE_STRIP);

		for(int i=0;i<ivx.size();i++)
		{
				glVertex2f(ivx[i],ivy[i]);
			
		}

			glEnd();
			glDisable(GL_LINE_STIPPLE);
			glColor3f(1.0f,1.0f,1.0f);

		//end guide lines
	}

	
	// show key with slide effect
	if(showKey)
		{
			if(showKeyPos<9)
			{
				showKeyPos++;
			}

			glColor3f(0.0f,0.5f,1.0f);
			renderBitmapCharacter(9+showKeyPos,17,0,GLUT_BITMAP_8_BY_13,"Quad/Cubic B-Spline Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,15.0f);
			glEnd();
			renderCircle(5+showKeyPos,15,2);

			glColor3f(1.0f,0.3f,0.3f);
			renderBitmapCharacter(9+showKeyPos,27,0,GLUT_BITMAP_8_BY_13,"Subdivision B-Spline Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,25.0f);
			glEnd();
			renderCircle(5+showKeyPos,25,2);


			glColor3f(1.0f,0.0f,1.0f);
			renderBitmapCharacter(9+showKeyPos,37,0,GLUT_BITMAP_8_BY_13,"Parametric Biezer Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,35.0f);
			glEnd();
			renderCircle(5+showKeyPos,35,2);

			glColor3f(0.6f,0.6f,1.0f);
			renderBitmapCharacter(9+showKeyPos,47,0,GLUT_BITMAP_8_BY_13,"Cubic Biezer Curve C0/C1 ");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,45.0f);
			glEnd();
			renderCircle(5+showKeyPos,45,2);

			glColor3f(1.0f,1.0f,0.3f);
			renderBitmapCharacter(9+showKeyPos,57,0,GLUT_BITMAP_8_BY_13,"Quad Biezer Curve C0/C1 ");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,55.0f);
			glEnd();
			renderCircle(5+showKeyPos,55,2);

			glColor3f(1.0f,0.5f,0.0f);
			renderBitmapCharacter(9+showKeyPos,67,0,GLUT_BITMAP_8_BY_13,"SubDivision Biezer Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,65.0f);
			glEnd();
			renderCircle(5+showKeyPos,65,2);
	}
	else
	{
	if(showKeyPos>-50)
			{
				showKeyPos--;
			}
			glColor3f(0.0f,0.5f,1.0f);
			renderBitmapCharacter(9+showKeyPos,17,0,GLUT_BITMAP_8_BY_13,"Quad/Cubic B-Spline Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,15.0f);
			glEnd();
			renderCircle(5+showKeyPos,15,2);

			glColor3f(1.0f,0.3f,0.3f);
			renderBitmapCharacter(9+showKeyPos,27,0,GLUT_BITMAP_8_BY_13,"Subdivision B-Spline Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,25.0f);
			glEnd();
			renderCircle(5+showKeyPos,25,2);


			glColor3f(1.0f,0.0f,1.0f);
			renderBitmapCharacter(9+showKeyPos,37,0,GLUT_BITMAP_8_BY_13,"Parametric Biezer Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,35.0f);
			glEnd();
			renderCircle(5+showKeyPos,35,2);

			glColor3f(0.6f,0.6f,1.0f);
			renderBitmapCharacter(9+showKeyPos,47,0,GLUT_BITMAP_8_BY_13,"Cubic Biezer Curve C0/C1 ");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,45.0f);
			glEnd();
			renderCircle(5+showKeyPos,45,2);

			glColor3f(1.0f,1.0f,0.3f);
			renderBitmapCharacter(9+showKeyPos,57,0,GLUT_BITMAP_8_BY_13,"Quad Biezer Curve C0/C1 ");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,55.0f);
			glEnd();
			renderCircle(5+showKeyPos,55,2);

			glColor3f(1.0f,0.5f,0.0f);
			renderBitmapCharacter(9+showKeyPos,67,0,GLUT_BITMAP_8_BY_13,"SubDivision Biezer Curve");
			glBegin(GL_POINTS);
			glVertex2f(5.0f+showKeyPos,65.0f);
			glEnd();
			renderCircle(5+showKeyPos,65,2);

	}

	glColor3f(1.0f,1.0f,1.0f);
	glPopMatrix();
	perspective_projection();


	//for realtime interaction
	if(autoUpdate)
	{
		genCurve();
	}

	

	glutSwapBuffers();

	glFlush();

}

void processMenuEvents(int option) {

	switch (option) {
		case SELECT : 
			{
				if(ivx.size()>0)
				{
					select_on = 1;	
					curve.clear();

					for(int i = 0; i < ivx.size(); i++)
					{
						if(	(abs( ivx[i]-x_cord ) < 4) &&(abs( ivy[i]-y_cord ) < 4))
						{
							point_id = i;
							
						}
					}
				}

			} break;
		case DUPLICATE : 
			{
				if(ivx.size()>0)
				{
					duplicate_on = 1	;

					for(int i = 0; i < ivx.size(); i++)
					{
						if(	(abs( ivx[i]-x_cord ) < 4) &&(abs( ivy[i]-y_cord ) < 4))
						{
							ivx.push_back(0);
							ivy.push_back(0);
							for(int j = ivx.size()-1; j>i+1; j--)
							{
								ivx[j]=ivx[j-1];
								ivy[j]=ivy[j-1];
							}
							ivx[i+1]=ivx[i];
							ivy[i+1]=ivy[i];


							point_id = i+1;
							break;
						}
					}
				}

			} break;
		case DELETE : 
			{

				for(int i = 0; i < ivx.size(); i++)
				{
					if(	(abs( ivx[i]-x_cord ) < 4) &&(abs( ivy[i]-y_cord ) < 4))
					{
						if( i == ivx.size()-1)
						{
							curve.clear();
							ivx.pop_back();
							ivy.pop_back();
							break;
						}
						else
						{
							for(int j = i; j<ivx.size()-1; j++)
							{
								ivx[j]=ivx[j+1];
								ivy[j]=ivy[j+1];
							}
							curve.clear();
							ivx.pop_back();
							ivy.pop_back();
							break;
						}
						
					}
				}

			} break;
	case DELETE_ALL : 
			{	
				curve.clear();
				ivx.clear();
				ivy.clear();
				display_3D =1;

			} break;

	case REFRESH : 
			{	
				genCurve();
			} break;
	}
}
void createMenus() {

	int submenu;

	submenu = glutCreateMenu(processMenuEvents);

	glutAddMenuEntry("Select",SELECT);
	glutAddMenuEntry("Duplicate",DUPLICATE);
	glutAddMenuEntry("Delete",DELETE);
	glutAddMenuEntry("Delete All",DELETE_ALL);
	glutAddMenuEntry("Refresh",REFRESH);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}


//Initializes 3D rendering
void initRendering() {
	//Makes 3D drawing work when something is in front of something else
	glEnable(GL_DEPTH_TEST);

}

//Called when the window is resized
void handleResize(int w, int h) {
	//Tell OpenGL how to convert from coordinates to pixel values
	glViewport(0, 0, w, h);
	
	x_pixel = w;
	y_pixel = h;	
	
		//motion blur
		glClearAccum(0.0, 0.0, 0.0, 1.0);
		glClear(GL_ACCUM_BUFFER_BIT);
		//motion blur

	glMatrixMode(GL_PROJECTION); //Switch to setting the camera perspective
	
	//Set the camera perspective
	glLoadIdentity(); //Reset the camera
	gluPerspective(45.0,                  //The camera angle
				   (double)w / (double)h, //The width-to-height ratio
				   0.0,                   //The near z clipping coordinate
				   200.0);                //The far z clipping coordinate
}

//Draws the 3D scene
void drawGlutIdle( void )
{
  /* According to the GLUT specification, the current window is 
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
  if ( glutGetWindow() != main_window ) 
    glutSetWindow(main_window);  

  glutPostRedisplay();	

  glui->sync_live();

}


void update(int value) {

	if(display_3D)
	{
		_angle += 2.0f;
		if (_angle > 360)
		{
			_angle -= 360;
		}
		if(sketchEffect)
		{		
				introID = 2;
				glClearColor(1.0f,1.0f,1.0f,1.0f);
				glEnable (GL_DEPTH_TEST); 
				GLfloat density = 5.7;
				GLfloat fogColor[4] = {0.0, 0.0, 0.0, 1.0}; 
				glEnable (GL_FOG);
				glFogi (GL_FOG_MODE, GL_EXP2); //set the fog mode to GL_EXP2
				glFogfv (GL_FOG_COLOR, fogColor); //set the fog color toour color chosen above
				glFogf (GL_FOG_DENSITY, density); //set the density to thevalue above
				glHint (GL_FOG_HINT, GL_NICEST); // set the fog to look thenicest, may slow down on older cards
		}
		else
		{		glClearColor(0.0f,0.0f,0.0f,1.0f);
				glDisable (GL_DEPTH_TEST);
				glDisable (GL_FOG);
		}
	}
	else
	{
		glClearColor(0.0f,0.0f,0.0f,1.0f);
		glDisable (GL_DEPTH_TEST);
		glDisable (GL_FOG);
	}

	tStrength=1/(double)strength;

	glutPostRedisplay(); //Tell GLUT that the display has changed
	
	//Tell GLUT to call update again in 25 milliseconds
	glutTimerFunc(25, update, 0);



}

//Called when a key is pressed
void handleKeypress(unsigned char key, int x, int y) {   //The key that was pressed n The current mouse coordinates
	switch (key) {
		case 27: //Escape key
			{
				exit(0); //Exit the program
				break;
			}
		case 13:						// enter(carraige return) key
			{	
				if(fullscreen == 0)
				{
					glutFullScreen();
					fullscreen = 1;
				}
				else
				{
					glutReshapeWindow(878,540);
					glutPositionWindow(200, 200);
					fullscreen = 0;
				}
					break;
			}
		default :
			{
				x_cord = x;
				y_cord = y;
				break;
			}
	}
}
void  handleMotion(int x,int y)
{
	x_cord = x;
	y_cord = y;

	if(select_on && point_id<ivx.size())
		{
				ivx[point_id] = x_cord;
				ivy[point_id] = y_cord;
		}
	if(duplicate_on)
		{

			int snapx,snapy;
			snapx = ivx[point_id-1] ;
			snapy = ivy[point_id-1];
			//code to snap close points
			if(abs(snapx-x)<5 && abs(snapy-y)<5)
			{
				ivx[point_id] =snapx;
				ivy[point_id]=snapy;
				snap = 1;
			}
			else
			{
				ivx[point_id] = x_cord;
				ivy[point_id] = y_cord;
				snap = 0;
			}

		}

}

void handleMouse(int button, int button_state, int x, int y )
{
  if ( button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN ) 
  {
	display_3D = 0;
    x_cord = x;
    y_cord = y;
	toggle = 0;

	if(select_on)
	{
		select_on = 0;
	}
	else if(duplicate_on)
		 {		
			duplicate_on = 0;
			snap = 0;
		}

	else
	{
	
		ivx.push_back(x);
		ivy.push_back(y);

	}
  }

  //genCurve();

}

int main(int argc, char** argv) {

	//Initialize GLUT
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(878, 490);
	glutInitWindowPosition(200, 200);

	initRendering(); //Initialize rendering
	
	main_window = glutCreateWindow( "GLUT OUT" );


  GLUI_Master.set_glutReshapeFunc( handleResize );  
  GLUI_Master.set_glutKeyboardFunc( handleKeypress );
  GLUI_Master.set_glutSpecialFunc( NULL );
  GLUI_Master.set_glutMouseFunc( handleMouse );


	glutReshapeWindow(878,540);
	glutPositionWindow(200, 200);
	glutDisplayFunc(drawScene);
	//glutReshapeFunc(handleResize);
	//glutKeyboardFunc(handleKeypress);
	//glutSpecialFunc( NULL );
	//glutMouseFunc(handleMouse);	
	glutMotionFunc(handleMotion);	
	glutPassiveMotionFunc(handleMotion); 


/*
	//-----------------------------------------------------------
	//GLUI PART
	
	//-----------------------------------------------------------
	*/
	//-----------------------------------------------------------
	//GLUI PART
	

  /****************************************/
  /*         Here's the GLUI code         */
  /****************************************/
glui = GLUI_Master.create_glui_subwindow( main_window,GLUI_SUBWINDOW_RIGHT );

  //glui = GLUI_Master.create_glui( "GLUI OUT", 0, 200, 200 ); /* name, flags, x, and y */
   char tmp[50];
  sprintf( tmp ,"GLUI version: %3.2f\n", GLUI_Master.get_version() );
  new GLUI_StaticText( glui, tmp ); 
  new GLUI_Separator( glui );

  new GLUI_StaticText( glui, "Lab 1 Curve Modelling" ); 
	hudOptions = new GLUI_Panel(glui, "Basic Menu and HUD options" );
		// Add the Draw Check box to the 'Object Properties' Panel

		glui->add_checkbox_to_panel (hudOptions, "Show intro", &display_3D );
		glui->add_checkbox_to_panel (hudOptions, "Sketch Effect", &sketchEffect );
		spinnerIntroID  = new GLUI_Spinner( hudOptions, "Intro ID:", &introID);
		spinnerIntroID->set_int_limits( 1, 3 );

		glui->add_checkbox_to_panel (hudOptions, "Show Key", &showKey );
		glui->add_checkbox_to_panel (hudOptions, "Show Guide Lines", &showGuide );


	curveOptions = new GLUI_Panel(glui,"Curve Type:");
			glui->add_checkbox_to_panel (curveOptions, "Auto Update Curve", &autoUpdate );
			glui->add_checkbox_to_panel (curveOptions,"Show Division Points", &showPoints );
			new GLUI_Separator( curveOptions );
			radio_group = new GLUI_RadioGroup(curveOptions,&RADIO);
 			glui->add_radiobutton_to_group( radio_group, "Quad B-Spline" );		
			glui->add_radiobutton_to_group( radio_group, "Cube B-Spline" );	
			glui->add_radiobutton_to_group( radio_group, "SubDivision B-spline" );
			glui->add_radiobutton_to_group( radio_group, "Parametric Bezier" );			
			glui->add_radiobutton_to_group( radio_group, "Cubic Bezier with C0 Continuity" );
			glui->add_radiobutton_to_group( radio_group, "Cubic Bezier with C1 Continuity" );
			glui->add_radiobutton_to_group( radio_group, "Quad Bezier with C0 Continuity" );
			glui->add_radiobutton_to_group( radio_group, "Quad Bezier with C1 Continuity" );
			glui->add_radiobutton_to_group( radio_group, "SubDivision Bezier" );
			new GLUI_Separator( curveOptions );
	spinner  = new GLUI_Spinner( curveOptions, "SubDivision Segments:", &segments);
	spinner->set_int_limits( 1, 20 );
	spinner2  = new GLUI_Spinner( curveOptions, "Parametric t[1,0] strength:", &strength);
	spinner2->set_int_limits( 1, 200 );



	cordDisplay=new GLUI_Panel(glui,"Cordinates");
	GLUI_EditText *x_cord_text = new GLUI_EditText( cordDisplay, "X:", &x_cord );
	x_cord_text->disable();
	GLUI_EditText *y_cord_text = new GLUI_EditText( cordDisplay, "Y:", &y_cord );
	y_cord_text->disable();

  //A 'quit' button //
  new GLUI_Button(glui, "Quit", 0,(GLUI_Update_CB)exit );

  // Link windows to GLUI, and register idle callback //
  glui->set_main_gfx_window( main_window );


 // We register the idle callback with GLUI, not with GLUT//
  GLUI_Master.set_glutIdleFunc( drawGlutIdle );

  // Regular GLUT main loop ///

	createMenus();	
	glutTimerFunc(25, update, 0);
	glutMainLoop(); //Start the main loop.  glutMainLoop doesn't return.

	//This line is never reached  
  	return 0; 

}