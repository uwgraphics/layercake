class Wedge {
void setup()
{
  size(600,600);
  stroke(0);
}

float wedgeAngle  = PI * .2;
float wedgeSize = 400;
private int wedgeX = 300;
private int wedgeY = 550;

float right = .7;
float left = .3;
float top = .7;
float bottom = .3;

int mygray=200;

//these are for the single control mouse point right now
int mx = wedgeX;
int my = wedgeY;
int mouseState = 0;

public Wedge(){
  }
  
public Wedge(int xpos,int ypos, int wsize){
    this.wedgeX = xpos; 
    this.mx = wedgeX;
    this.wedgeY = ypos; 
    this.my = wedgeY;
    this.wedgeSize = wsize;
  }  

public void setWedgeX(int xpos){
    this.wedgeX = xpos; 
    this.mx = wedgeX;
  }
  
public void setWedgeY(int ypos){
    this.wedgeY = ypos; 
    this.my = wedgeY;
  }

float lerpAngle(float pct)
{
  return (-PI/2- (1-2*pct) * wedgeAngle);
}

void myarc(float p1,float p2,float r)
{
  arc(wedgeX,wedgeY,wedgeSize*r,wedgeSize*r,lerpAngle(p1),lerpAngle(p2));
}

//these control the 2 yellow button positions
float posX(float th, float r)
{
  return wedgeX + r*wedgeSize*cos(lerpAngle(th))/2;
}
float posY(float th, float r)
{
  return wedgeY + r*wedgeSize*sin(lerpAngle(th))/2;
}

void draw()
{
  pushStyle();
  pushMatrix();
  noStroke();
  colorMode(HSB,1.0,1.0,1.0);
  
  ellipseMode(CENTER);
    fill(1, 1, .76);
  myarc(right,1,1);
  fill(1, .84, .76);
  myarc(right,1,top);
  fill(1, 0.67, .76);
  myarc(right,1,top*4/5+bottom/5);
  fill(1, 0.50, .77);
  myarc(right,1,top*3/5+bottom*2/5);
  fill(1, 0.34, .77);
  myarc(right,1,top*2/5+bottom*3/5);
  fill(1, 0.17, .78);
  myarc(right,1,top/5  +bottom*4/5);
  
  //right 90 
  fill(.97,.95,.6);
  myarc(left/10 + right*9/10,right,1);
  fill(.97,.79,.63);
  myarc(left/10 + right*9/10,right,top);
  fill(.97,.64,.66);//.2,.6,.7);
  myarc(left/10 + right*9/10,right,top*4/5+bottom/5);
  fill(.97,.48,.70);//.2,.6,.7);
  myarc(left/10 + right*9/10,right,top*3/5+bottom*2/5);
  fill(.97,.32,.74);//.2,.6,.7);
  myarc(left/10 + right*9/10,right,top*2/5+bottom*3/5);
  fill(.97,.16,.78);//.2,.6,.7);
  myarc(left/10 + right*9/10,right,top/5+bottom*4/5);
  
  //right 80 
  fill(.99,.90,.7);
  myarc(left*2/10 + right*8/10,left/10 + right*9/10,1);
  fill(.99,.75,.71);
  myarc(left*2/10 + right*8/10,left/10 + right*9/10,top);
  fill(.99,.60,.73);//.2,.6,.7);
  myarc(left*2/10 + right*8/10,left/10 + right*9/10,top*4/5+bottom/5);
  fill(.99,.45,.74);//.2,.6,.7);
  myarc(left*2/10 + right*8/10,left/10 + right*9/10,top*3/5+bottom*2/5);
  fill(.99,.30,.76);//.2,.6,.7);
  myarc(left*2/10 + right*8/10,left/10 + right*9/10,top*2/5+bottom*3/5);
  fill(.99,.15,.78);//.2,.6,.7);
  myarc(left*2/10 + right*8/10,left/10 + right*9/10,top/5+bottom*4/5);
  
  //right 70
  fill(.02,.86,0.79);
  myarc(left*3/10 + right*7/10,        left*2/10 + right*8/10,  1);
  fill(.02,.72,0.73);
  myarc(left*3/10 + right*7/10,        left*2/10 + right*8/10,  top);
  fill(.02,.58,0.74);//.2,.6,.7);
  myarc(left*3/10 + right*7/10,        left*2/10 + right*8/10,  top*4/5+bottom/5);
  fill(.02,.44,0.75);//.2,.6,.7);
  myarc(left*3/10 + right*7/10,        left*2/10 + right*8/10,  top*3/5+bottom*2/5);
  fill(.02,.29,0.76);//.2,.6,.7);
  myarc(left*3/10 + right*7/10,        left*2/10 + right*8/10,  top*2/5+bottom*3/5);
  fill(.02,.14,0.78);//.2,.6,.7);
  myarc(left*3/10 + right*7/10,        left*2/10 + right*8/10,  top/5+bottom*4/5);
  
  //right 60
  fill(.05,.81,0.89);
  myarc(left*4/10 + right*6/10,        left*3/10 + right*7/10,  1);
  fill(.05,.67,0.86);
  myarc(left*4/10 + right*6/10,        left*3/10 + right*7/10,  top);
  fill(.05,.54,0.84);//.2,.6,.7);
  myarc(left*4/10 + right*6/10,        left*3/10 + right*7/10,  top*4/5+bottom/5);
  fill(.05,.40,0.82);//.2,.6,.7);
  myarc(left*4/10 + right*6/10,        left*3/10 + right*7/10,  top*3/5+bottom*2/5);
  fill(.05,.27,0.80);//.2,.6,.7);
  myarc(left*4/10 + right*6/10,        left*3/10 + right*7/10,  top*2/5+bottom*3/5);
  fill(.05,.13,0.78);//.2,.6,.7);
  myarc(left*4/10 + right*6/10,        left*3/10 + right*7/10,  top/5+bottom*4/5);
  
  //mid
  fill(.069, .76,0.89);
  myarc(left*5/10 + right*5/10,        left*4/10 + right*6/10,  1);
  fill(.069 , .63,0.86);
  myarc(left*5/10 + right*5/10,        left*4/10 + right*6/10,  top);
  fill(.069 , .50,0.84);//.2,.6,.7);
  myarc(left*5/10 + right*5/10,        left*4/10 + right*6/10,  top*4/5+bottom/5);
  fill(.069 , .37,0.82);//.2,.6,.7);
  myarc(left*5/10 + right*5/10,        left*4/10 + right*6/10,  top*3/5+bottom*2/5);
  fill(.069 , .25,0.80);//.2,.6,.7);
  myarc(left*5/10 + right*5/10,        left*4/10 + right*6/10,  top*2/5+bottom*3/5);
  fill(.069 , .12,0.78);//.2,.6,.7);
  myarc(left*5/10 + right*5/10,        left*4/10 + right*6/10,  top/5+bottom*4/5);
  
  //left 40 
  fill(.09,.65,.99);
  myarc(left*6/10 + right*4/10,        left*5/10 + right*5/10,  1);
  fill(.09,.54,.95);
  myarc(left*6/10 + right*4/10,        left*5/10 + right*5/10,  top);
  fill(.09,.43,.90);//.2,.6,.7);
  myarc(left*6/10 + right*4/10,        left*5/10 + right*5/10,  top*4/5+bottom/5);
  fill(.09,.32,.86);//.2,.6,.7);
  myarc(left*6/10 + right*4/10,        left*5/10 + right*5/10,  top*3/5+bottom*2/5);
  fill(.09,.22,.82);//.2,.6,.7);
  myarc(left*6/10 + right*4/10,        left*5/10 + right*5/10,  top*2/5+bottom*3/5);
  fill(.09,.11,.78);//.2,.6,.7);
  myarc(left*6/10 + right*4/10,        left*5/10 + right*5/10,  top/5+bottom*4/5);
    
  //left 30 
  fill(.108,.54,.99);
  myarc(left*7/10 + right*3/10,        left*6/10 + right*4/10,  1);
  fill(.108,.45,.95);
  myarc(left*7/10 + right*3/10,        left*6/10 + right*4/10,  top);
  fill(.108,.36,.90);//.2,.6,.7);
  myarc(left*7/10 + right*3/10,        left*6/10 + right*4/10,  top*4/5+bottom/5);
  fill(.108,.27,.86);//.2,.6,.7);
  myarc(left*7/10 + right*3/10,        left*6/10 + right*4/10,  top*3/5+bottom*2/5);
  fill(.108,.18,.82);//.2,.6,.7);
  myarc(left*7/10 + right*3/10,        left*6/10 + right*4/10,  top*2/5+bottom*3/5);
  fill(.108,.09,.78);//.2,.6,.7);
  myarc(left*7/10 + right*3/10,        left*6/10 + right*4/10,  top/5+bottom*4/5);
  
  //left 20  
  fill(.128,.42,1);
  myarc(left*8/10 + right*2/10,        left*7/10 + right*3/10,  1);
  fill(.128,.35,.96);
  myarc(left*8/10 + right*2/10,        left*7/10 + right*3/10,  top);
  fill(.128,.28,.91);//.2,.6,.7);
  myarc(left*8/10 + right*2/10,        left*7/10 + right*3/10,  top*4/5+bottom/5);
  fill(.128,.21,.87);//.2,.6,.7);
  myarc(left*8/10 + right*2/10,        left*7/10 + right*3/10,  top*3/5+bottom*2/5);
  fill(.128,.14,.82);//.2,.6,.7);
  myarc(left*8/10 + right*2/10,        left*7/10 + right*3/10,  top*2/5+bottom*3/5);
  fill(.128,.07,.78);//.2,.6,.7);
  myarc(left*8/10 + right*2/10,        left*7/10 + right*3/10,  top/5+bottom*4/5);
  
  //left 10  
  fill(.147,.31,1);
  myarc(left*9/10 + right/10,        left*8/10 + right*2/10,  1);
  fill(.147,.26,.96);
  myarc(left*9/10 + right/10,        left*8/10 + right*2/10,  top);
  fill(.147,.21,.91);//.2,.6,.7);
  myarc(left*9/10 + right/10,        left*8/10 + right*2/10,  top*4/5+bottom/5);
  fill(.147,.16,.87);//.2,.6,.7);
  myarc(left*9/10 + right/10,        left*8/10 + right*2/10,  top*3/5+bottom*2/5);
  fill(.147,.11,.82);//.2,.6,.7);
  myarc(left*9/10 + right/10,        left*8/10 + right*2/10,  top*2/5+bottom*3/5);
  fill(.147,.05,.78);//.2,.6,.7);
  myarc(left*9/10 + right/10,        left*8/10 + right*2/10,  top/5+bottom*4/5);  
  
  //left 
  fill(.167,.2,1);
  myarc(left,  left*9/10 + right/10,  1);
  fill(.167,.17,0.96);
  myarc(left,  left*9/10 + right/10,  top);
  fill(.167,.14,0.91);//.2,.6,.7);
  myarc(left,  left*9/10 + right/10,  top*4/5+bottom/5);
  fill(.167,.10,0.87);//.2,.6,.7);
  myarc(left,  left*9/10 + right/10,  top*3/5+bottom*2/5);
  fill(.167,.07,0.82);//.2,.6,.7);
  myarc(left,  left*9/10 + right/10,  top*2/5+bottom*3/5);
  fill(.167,.03,0.78);//.2,.6,.7);
  myarc(left,  left*9/10 + right/10,  top/5+bottom*4/5);
  
  fill(0,0,.7);
  myarc(left,1,bottom);
  
  colorMode(RGB,1.0,1.0,1.0);
  fill(1,1,0);
  int rad1 = 10;
  int rad2 = 10;
  if(mouseState==1){
    rad1 = 15;
  }
  else if(mouseState==2){
    rad2 = 15;
  }
  ellipse(posX(right,top),posY(right,top),rad1,rad1);
  ellipse(posX(left,bottom),posY(left,bottom),rad2,rad2);
  if(mouseState>0){
   // ellipse(mx,my,9,9);
  }
  popMatrix();
  popStyle();
}

void mousePressedAlt(float mX,float mY)
{
  if ((abs(posX(right,top)-mX) < 20) &&
      (abs(posY(right,top)-mY) < 20)) {
    mouseState = 1;
      }
  if ((abs(posX(left,bottom)-mX) < 20) &&
      (abs(posY(left,bottom)-mY) < 20)) {
    mouseState = 2;
      }
}
void mouseReleasedAlt()
{
  mouseState = 0;
}

void mouseDraggedAlt(float nx,float ny)
{
  mx = int(nx);
  my = int(ny);  

  float fx = float(mx-wedgeX);
  float fy = float(wedgeY-my);
  
  float r = sqrt(fx*fx+fy*fy)/wedgeSize*2;
  float a = atan2(fy,fx);
  
  a = 1-(1+(a-PI/2)/wedgeAngle)/2;
  
  if (a<0) a = 0;
  if (a>1) a = 1;
  if (r<0) r = 0;
  if (r>1) r = 1;
  
  if (mouseState == 1) {
    right = a;
    top = r;
    
    if (left>right) left=right;
    if (bottom>top) bottom=top;
  }
  if (mouseState == 2) {
    left = a;
    bottom = r;
    
    if(right<left) right=left;
    if(top<bottom) top=bottom;
  }
}


}

