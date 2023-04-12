

#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include<bits/stdc++.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<iostream>
#include <vector>
#include <string>
#include <sstream>
#include "1705006_classes.h"
#include "bitmap_image.hpp"
#define pi (2*acos(0.0))
#define PI (2*acos(0.0))

#define windowWidth 500
#define windowHeight 500
using namespace std;
string PATH="D:\\StudyTools\\L4T1\\Computer Graphics\\OPEN GL projects\\offline3_demo\\";

int grid_difference = 15;
int grid_len = 185;

//int wheel_radius = 35;
//int wheel_height = 18;
int wheel_segments = 40;
int angle_step = 2;
int distance_step = 2;
int count_bmp = 0;

//int draw_grid;
//int draw_axis;

int number_of_objs,number_of_pointligths,number_of_spotlights;

int drawgrid;
int drawaxes;


int imageWidth;
int imageHeight;
///extern VARS
vector<Object*> objects;
vector<PointLight*> point_lights;
vector<SpotLight*> spot_lights;
int recursion_level;



double DegreeToRad(double degree)
{
    return (degree*PI)/180;
}

///--------------------------FULLY CONTROLLABLE CAMERA--------------------///
class Camera
{
private:
    point position = point(120, 120, 40);
    point look = point(-1 / sqrt(2), -1 / sqrt(2), 0);
    point right = point(-1 / sqrt(2), 1 / sqrt(2), 0);
    point up = point(0, 0, 1);
public:

    ///getters
    const point& get_position() const {
		return position;
	}

	const point& get_look() const {
		return look;
	}
	const point& get_right() const {
		return right;
	}

	const point& get_up() const {
		return up;
	}
	///IN LOOK ==> ONLY UP & LOOK VECTOR CHANGES
    void Lookup()
    {
        ///change up and look direction
        point temp = up;
        up = up * cos(DegreeToRad(2)) - look * sin(DegreeToRad(2));
		look = look * cos(DegreeToRad(2)) + temp * sin(DegreeToRad(2));
    }
    void Lookdown() {
        ///change up and look direction
        point temp = up;
		up = up * cos(-DegreeToRad(2)) - look * sin(-DegreeToRad(2));
		look = look * cos(-DegreeToRad(2)) + temp * sin(-DegreeToRad(2));
	}

	///in look == ONLY RIGHT AND LOOK VECTORS CHANGES
	void Lookleft()
	{
	    ///change right and look direction
	    point temp = look;
	    look = look * cos(DegreeToRad(2)) - right * sin(DegreeToRad(2));
		right = right * cos(DegreeToRad(2)) + temp * sin(DegreeToRad(2));
	}
	void Lookright()
	{
	    ///change right and look direction
	    point temp = look;
	    look = look * cos(-DegreeToRad(2)) - right * sin(-DegreeToRad(2));
		right = right * cos(-DegreeToRad(2)) + temp * sin(-DegreeToRad(2));

	}

	///in tilt only RIGHT & UP vector will change LOOK vector won't change
	void Rotate_clockwise()
	{
	    ///change up and right direction
	    point temp = right;
        right = right * cos(DegreeToRad(2)) - up * sin(DegreeToRad(2));
		up = up * cos(DegreeToRad(2)) + temp * sin(DegreeToRad(2));
	}
	void Rotate_anticlockwise()
	{
	    ///change up and right direction
	    point temp = right;
	    right = right * cos(-DegreeToRad(2)) - up * sin(-DegreeToRad(2));
		up = up * cos(-DegreeToRad(2)) + temp * sin(-DegreeToRad(2));

	}
	void Moveup()
	{
	    position += 2* up;

	}
	void Movedown()
	{
	    position -= 2 * up;

	}
	void Moveleft()
	{
        position -= 2 * right;
	}
	void Moveright()
	{
        position += 2 * right;
	}
	void Moveforward() {
		position += 2 * look;
	}

	void Movebackward() {
		position -= 2 * look;
	}

};
Camera camera;  ///creating global camera object

///--------------------------FULLY CONTROLLABLE CAMERA--------------------///



///---------------------draw grid----------------///
void drawGrid() {
	if (drawgrid == 1) {
        glColor3f(0.7, 0.7, 0.7);
        auto len = (grid_len / grid_difference);
		glBegin(GL_LINES);
		{
			for (auto i = -len; i <= len; i++) {

				glVertex3f(i * grid_difference, -grid_len, 0);
				glVertex3f(i * grid_difference, grid_len, 0);


				glVertex3f(-grid_len, i * grid_difference, 0);
				glVertex3f(grid_len, i * grid_difference, 0);
			}
		}
		glEnd();
	}
}
///----------------draw axes--------------------///
void drawAxes() {
	if (drawaxes == 1) {
		glBegin(GL_LINES);
		{
			/// x-axis
			glColor3f(1, 0, 0);
			glVertex3f(200, 0, 0);
			glVertex3f(-200, 0, 0);

			/// y-axis
			glColor3f(0, 1, 0);
			glVertex3f(0, -200, 0);
			glVertex3f(0, 200, 0);

			/// z-axis
			glColor3f(0, 0, 1);
			glVertex3f(0, 0, 200);
			glVertex3f(0, 0, -200);
		}
		glEnd();
	}
}

std::string toString(auto &i) {
    std::stringstream ss;
    ss << i;

    return ss.str();
}



void capture()
{

    cout<<"------------------CAPTURE STARTS----------------"<<endl;
    bitmap_image image(imageWidth,imageHeight);
    for(int i=0;i<imageWidth;i++){
        for(int j=0;j<imageHeight;j++){
            image.set_pixel(i,j,0,0,0);
        }
    }

    string temp_string = "1705006_out" + toString(count_bmp)+".bmp";

    string img_path=PATH+temp_string;
    image.save_image(img_path);;
    double planeDistance=(windowHeight/2)/(tan(80*(pi/360)));


    double tempX = camera.get_position().x + (camera.get_look().x * planeDistance) - camera.get_right().x * (windowWidth/2.0) + camera.get_up().x * (windowHeight/2.0);
    double tempY = camera.get_position().y + (camera.get_look().y * planeDistance) - camera.get_right().y * (windowWidth/2.0) + camera.get_up().y * (windowHeight/2.0);
    double tempZ = camera.get_position().z + (camera.get_look().z * planeDistance) - camera.get_right().z * (windowWidth/2.0) + camera.get_up().z * (windowHeight/2.0);

    struct point topleft(tempX,tempY,tempZ);



    double du=(windowWidth*1.0)/imageWidth;
    double dv=(windowHeight*1.0)/imageHeight;


    topleft.x = topleft.x + (camera.get_right().x*(0.5*du)) - (camera.get_up().x*0.5*dv);
    topleft.y = topleft.y + (camera.get_right().y*(0.5*du)) - (camera.get_up().y*0.5*dv);
    topleft.z = topleft.z + (camera.get_right().z*(0.5*du)) - (camera.get_up().z*0.5*dv);

    int nearest=-1;
    double temp_t;

    for(int i=0;i<imageWidth;i++){
        for(int j=0;j<imageHeight;j++){
            struct point curPixel;
            curPixel.x = topleft.x + (camera.get_right().x * (i * du)) - (camera.get_up().x * j * dv);
            curPixel.y = topleft.y + (camera.get_right().y * (i * du)) - (camera.get_up().y * j * dv);
            curPixel.z = topleft.z + (camera.get_right().z * (i * du)) - (camera.get_up().z * j * dv);
            struct point RayCurlpixelEye(curPixel.x - camera.get_position().x,curPixel.y - camera.get_position().y,curPixel.z - camera.get_position().z);
            RayCurlpixelEye = unit_vector(RayCurlpixelEye);
            Ray CastRay(camera.get_position(),RayCurlpixelEye);


            double tMin=999999;
            double *color=new double[3];


            double *dummyColor=new double[3];
            color[0]=color[1]=color[2]=0.0;
            for(int k=0;k<objects.size();k++){
                temp_t=objects[k]->intersect(&CastRay,dummyColor,0);

                if(temp_t<=0){
                    continue;
                }
                else if(temp_t<tMin){

                    ///update tMIN & NEAREST OBJ
                    tMin=temp_t;
                    nearest=k;


                }
            }
            ///IF THERE EXISTS A NEAREST OBJ
            if(nearest!=-1){
                tMin=objects[nearest]->intersect(&CastRay,color,1);
                image.set_pixel(i,j,color[0]*255,color[1]*255,color[2]*255);

            }


            delete[] dummyColor;
            delete[] color;
        }

    }

    image.save_image(img_path);;
    cout<<"------------------CAPTURE ENDS----------------"<<endl;
    count_bmp++;



}

void loadData()
{
    ifstream in;

    string input_path=PATH+"scene.txt";
    in.open(input_path);

    if(!in.is_open()) {
        perror("Error!!!!file can't be opened");
        exit(EXIT_FAILURE);
    }
    in>>recursion_level;
    in>>imageHeight;
    imageWidth=imageHeight;
    in>>number_of_objs;
    Object *temp;
    for(int i=0;i<number_of_objs;i++){
        string readLine;
        in>>readLine;
        if(readLine.compare("sphere")==0|| readLine.compare("SPHERE")==0 || readLine.compare("Sphere")==0){
            double x,y,z,radius;
            in>>x>>y>>z;
            in>>radius;
            struct point center(x,y,z);

            temp=new Sphere(center,radius);
            double r,g,b;
            in>>r>>g>>b;
            temp->setColor(r,g,b);

            double AMB,DIF,SPEC,REF_COEFF;
            in>>AMB>>DIF>>SPEC>>REF_COEFF;

            temp->setCoEfficients(AMB,DIF,SPEC,REF_COEFF);
            int shine;
            in>>shine;
            temp->setShine(shine);
            objects.push_back(temp);
        }
        else if(readLine.compare("triangle")==0 || readLine.compare("TRIANGLE")==0 || readLine.compare("Triangle")==0){
            struct point x1;
            struct point y1;
            struct point z1;
            in>>x1.x>>x1.y>>x1.z;
            in>>y1.x>>y1.y>>y1.z;
            in>>z1.x>>z1.y>>z1.z;

            temp=new Triangle(x1,y1,z1);
            double r,g,b;
            in>>r>>g>>b;
            temp->setColor(r,g,b);
            double AMB,DIF,SPEC,REF_COEFF;
            in>>AMB>>DIF>>SPEC>>REF_COEFF;

            temp->setCoEfficients(AMB,DIF,SPEC,REF_COEFF);
            int shine;
            in>>shine;
            temp->setShine(shine);
            objects.push_back(temp);
        }
        else if(readLine.compare("general")==0 || readLine.compare("GENERAL")==0 || readLine.compare("General")==0){
            double A,B,C,D,E,F,G,H,I,J;


            in>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            double ref_x,ref_y,ref_z,length,width,height;
            in>>ref_x>>ref_y>>ref_z>>length>>width>>height;
            struct point referencePoint(ref_x,ref_y,ref_z);

            temp=new GeneralQuadratic(A,B,C,D,E,F,G,H,I,J,length,width,height,referencePoint);
            double r,g,b;
            in>>r>>g>>b;
            temp->setColor(r,g,b);
            double AMB,DIF,SPEC,REF_COEFF;
            in>>AMB>>DIF>>SPEC>>REF_COEFF;

            temp->setCoEfficients(AMB,DIF,SPEC,REF_COEFF);
            int shine;
            in>>shine;
            temp->setShine(shine);
            objects.push_back(temp);
        }
    }
    in>>number_of_pointligths;
    PointLight *tmp;
    for(int j=0;j<number_of_pointligths;j++){

        double p_x,p_y,p_z;
        in>>p_x>>p_y>>p_z;

        struct point point_light_position(p_x,p_y,p_z);

        tmp=new PointLight(point_light_position);
        double r,g,b;
        in>>r>>g>>b;
        tmp->setColor(r,g,b);
        point_lights.push_back(tmp);
    }
    in>>number_of_spotlights;
    SpotLight *tmp2;
    for(int j = 0;j < number_of_spotlights;j++)
    {
        PointLight *pl;
        double p_x,p_y,p_z;
        in>>p_x>>p_y>>p_z;

        struct point point_light_position(p_x,p_y,p_z);

        pl=new PointLight(point_light_position);
        double r,g,b;
        in>>r>>g>>b;
        pl->setColor(r,g,b);
        double dirX,dirY,dirZ;
        in>>dirX>>dirY>>dirZ;
        struct point light_dir(dirX,dirY,dirZ);
        //tmp2->light_dir = light_dir;
        //tmp2->point_light = pl;
        double cutoff_angle;
        in>>cutoff_angle;
        tmp2 = new SpotLight(pl,light_dir,cutoff_angle);
        //tmp2->cutoff_angle = cutoff_angle;
        spot_lights.push_back(tmp2);
    }
    temp=new Floor(1000,20);
    temp->setCoEfficients(0.34,0.25,0.2,0.2);
    temp->setShine(5);
    objects.push_back(temp);
    ///---------------------implement spotlights here--------------///
    ///in>>number_of_spotlights;


}
void Free_Memory()
{
    for(int i=0;i<objects.size();i++){
        delete objects[i];
    }
    for(int i=0;i<point_lights.size();i++){
        delete point_lights[i];
    }
    for(int i = 0;i< spot_lights.size();i++)
    {
        delete spot_lights[i];
    }
    objects.clear();
    point_lights.clear();
    spot_lights.clear();
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();
	///LOOK(CAMERA_POS  --        LOOKUP_POS   --       UP_POS)
	gluLookAt(camera.get_position().x, camera.get_position().y, camera.get_position().z, camera.get_position().x + camera.get_look().x,
		camera.get_position().y + camera.get_look().y, camera.get_position().z + camera.get_look().z, camera.get_up().x,
		camera.get_up().y, camera.get_up().z);



	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	//drawGrid();
    for(int i=0;i<objects.size();i++){
        objects[i]->draw();
    }
    for(int i=0;i<point_lights.size();i++){
        point_lights[i]->draw();
    }
    for(int i = 0;i<spot_lights.size();i++)
    {
        spot_lights[i]->point_light.draw();
    }


	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){


	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;




	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}





void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP

			}
			else if(state==GLUT_UP){

			}
			break;

		case GLUT_RIGHT_BUTTON:
		    if(state == GLUT_DOWN){
                drawaxes=1-drawaxes;
		    }

			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}


void keyboardListener(unsigned char key, int x, int y) {
	switch (key) {

	case '1':
		camera.Lookleft();
		break;
	case '2':
		camera.Lookright();
		break;

	case '3':
		camera.Lookup();
		break;
	case '4':
		camera.Lookdown();
		break;

	case '5':
		camera.Rotate_clockwise();
		break;
	case '6':
		camera.Rotate_anticlockwise();
		break;
    case '0':
        capture();
        break;

	default:
		break;
	}
}
void specialKeyListener(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_DOWN:        //down arrow key
		camera.Movebackward();
		break;
	case GLUT_KEY_UP:        // up arrow key
		camera.Moveforward();
		break;

	case GLUT_KEY_RIGHT:
		camera.Moveright();
		break;
	case GLUT_KEY_LEFT:
		camera.Moveleft();
		break;

	case GLUT_KEY_PAGE_UP:
		camera.Moveup();
		break;
	case GLUT_KEY_PAGE_DOWN:
		camera.Movedown();
		break;




	default:
		break;
	}
}





/* Program entry point */

int main(int argc, char *argv[])
{
    loadData();
    atexit(Free_Memory);
	glutInit(&argc,argv);
	glutInitWindowSize(windowWidth,windowHeight);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("--------------------------RAY TRACING---------------------------");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();
	 	//The main loop of OpenGL

	return 0;


}
