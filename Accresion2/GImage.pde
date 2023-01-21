class Pigment {
  int alpha;
  int red;
  int green;
  int blue;
  Pigment() {
  }
  Pigment(int r, int g, int b) {
    alpha=255;
    red=r;
    green=g;
    blue=b;
  }
  Pigment(int r, int g, int b,int a) {
    alpha=a;
    red=r;
    green=g;
    blue=b;
  }
}
class GImage {
  PImage bits;
  GImage(int w, int h) {
    bits=new PImage(w, h);
    Clear();
  }
  void Clear() {
    int w=bits.width*bits.height;
    for (int i=0; i<w; i++) {
      bits.pixels[i]=color(0, 0, 0, 0);
    }
  }
  void Clear(Pigment c) {

    int w=bits.width*bits.height;
    for (int i=0; i<w; i++) {
      bits.pixels[i]=color(c.red, c.green, c.blue, c.alpha);
    }
  }
  void ClearAdd(Pigment p, float light,float sat) {
    float r=p.red;
    float g=p.green;
    float b=p.blue;
    
    int w=bits.width*bits.height;
    float d=sat/(sat+light);
    for (int i=0; i<w; i++) {
      color c=bits.pixels[i];
      if(red(c)>0 ||green(c)>0 ||blue(c)>0 ){
      r=p.red;
      g=p.green;
      b=p.blue;
      r+=red(c)*light;
      g+=green(c)*light;
      b+=blue(c)*light;
      r*=d;
      g*=d;
      b*=d;
      bits.pixels[i]=color(r, g, b,p.alpha);
      }
    }
  }
  int At(int x, int y) {
    return Math.abs(x+y*bits.width)%(bits.width*bits.height);
  }
  int AtLimit(int x, int y) {
    if (x<0 || x>=bits.width || y<0 || y>=bits.height)return -1;
    return Math.abs(x+y*bits.width);
  }
  int Width() {
    return bits.width;
  }
  int Height() {
    return bits.height;
  }
  void AddAlphaPixel(int x, int y, int a, int r, int g, int b) {
    int at=At(x, y);
    color c=bits.pixels[at];
    r+=red(c);
    g+=green(c);
    b+=blue(c);
    r>>=1;
    g>>=1;
    b>>=1;
    bits.pixels[at]=color(r, g, b, 255);
  }
  void SetPixel(int x, int y, int a, int r, int g, int b) {
    int at=At(x, y);
    bits.pixels[at]=color(r, g, b);
  }
  void AddPixel(int x, int y, int a, int r, int g, int b) {
    int at=At(x, y);
    color c=bits.pixels[at];
    r+=red(c);
    g+=green(c);
    b+=blue(c);
    a+=alpha(c);
    r>>=1;
    g>>=1;
    b>>=1;
    a>>=1;
    bits.pixels[at]=color(r, g, b, a);
  }
  void SetPixel(float x, float y, int a, int r, int g, int b) {
    int at=AtLimit((int)x, (int)y);
    if (at>=0) {
      bits.pixels[at]=color(r, g, b);
    }
  }
  void AddPixel(float x, float y, int a, int r, int g, int b) {
    int at=AtLimit((int)x, (int)y);
    if (at>=0) {
      color c=bits.pixels[at];
      r+=red(c);
      g+=green(c);
      b+=blue(c);
      a+=alpha(c);
      r>>=1;
      g>>=1;
      b>>=1;
      a>>=1;
      bits.pixels[at]=color(r, g, b, a);
    }
  }
  void Update() {
    bits.updatePixels();
  }
  void Draw() {
    image(bits, 0, 0);
  }
  void Draw3D() {
    image(bits, -bits.width*0.5, -bits.height*0.5);
  }
}
