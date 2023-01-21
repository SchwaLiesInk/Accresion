/*
double f2(double x,double y,double a,double b){
  return a-x*x+b*y;
}
double sr=0.8+arcRandom(0.2);
double aa=1.4142135623730950488016887242097-arcRandom(0.03)+arcRandom(0.03);
double bb=1.6180339887498948482045868343656-aa;
*/
int counter=0;
float lys=7.0;
int epixel=250000;
//double au=1.5e11/6e6;//pixels
double scale=1.0/(lys*epixel);
double year=(24*60*60*356);
double day=(24*60*60);
double week=day*7;
double month=day*31;
double time=day;
double lightYear=(3.0e8)*year;
StarSystem sys;
float lspd=0.6;
//

double radius=3;
float spin=0;
int mx;  
int my;  
int trak=64;
StarSystemScreen ssScreen;
Vector coord;
void setup(){
size(800,800,P3D);
coord=new Vector();
coord.x=0;
coord.y=0;
coord.z=1000;
double MAXX=1000;
double MAXY=1000;
double MAXZ=1000;
arcRandomSeed((long)(MAXX+coord.x+(MAXY+coord.y)*MAXY+(MAXZ+coord.z)*MAXX*MAXY));
int n0=1+(int)(arcRandom(3)+arcRandom(2));//6;
int n1=n0+(int)(arcRandom(3)+arcRandom(2));//12;
    background(0);
    sys=new StarSystem(n1,10,null);
    ssScreen=new StarSystemScreen();
  textFont(createFont("Eurostile", 12));
}
void mousePressed() {
  ssScreen.mPressed();
}
void mouseReleased() {
  ssScreen.mReleased();
}
void mouseDragged() {
    ssScreen.mDragged();

}
void keyPressed() {

    ssScreen.kPressed();

}
void keyReleased() {

    ssScreen.kReleased();

}
void draw(){
  ssScreen.Keys();
  background(0);
  lights();
  ssScreen.Draw();
}
