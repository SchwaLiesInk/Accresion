import java.util.Random;
double uG=6.5e-11;
Random aRandom=new Random();
void arcRandomSeed(long s){
aRandom=new Random(s);
}
float arcRandom(float d){
return aRandom.nextFloat()*d;
}
double CalculateOrbitalVelocityMagnitude(double centerOfGravity, double dis)
{
  double ga = (uG * centerOfGravity) / (dis);
  return Math.sqrt(ga);
}
Vector CentrefugalForce(Vector velocity, Vector directionNormalToCenter, double distance)
{
  double v2 = velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z;
  double d2 = distance;
  double ratio = (v2 / d2);
  return new Vector(directionNormalToCenter.x * ratio, directionNormalToCenter.y * ratio, directionNormalToCenter.z * ratio);
}
double CentrefugalForceMagnitude(Vector velocity, double distance)
{
  double v2 = velocity.x * velocity.x + velocity.y * velocity.y + velocity.z * velocity.z;
  double d2 = distance;
  double ratio = (v2 / d2);
  return ratio;
}
double CentrefugalMagnitude(double velocity, double distance)
{
  double v2 = velocity * velocity;
  double d2 = distance;
  double ratio = (v2 / d2);
  return ratio;
}
double CircularVelocityMagnitude(double centeralMass, double distance)
{
  return Math.sqrt(uG * centeralMass / distance);
}
double EscapeVelocityMagnitude(double centeralMass, double distance)
{
  return Math.sqrt(2.0 * uG * centeralMass / distance);
}
double SurfaceGravityMagnitude(double centeralMass, double radius)
{
  return (uG *centeralMass) / (radius * radius);
}
double GravitationalForceMagnitude(double centeralMass1, double centeralMass2, double distance)
{
  return (uG * centeralMass1 * centeralMass2) / (distance * distance);
}
double MomentOfInertia(double mass, double distance)
{
  return mass * (distance * distance);
}
double MomentOfInertialSpin(double mass, double radius)
{
  return 0.4 * mass * (radius * radius);
}
double CurvedKineticEnergy(double momentOfIneria, double angularVelocity)
{
  return momentOfIneria * angularVelocity * angularVelocity * 0.5;
}
double KineticEnergy(double mass, double velocity)
{
  return mass * velocity * velocity * 0.5;
}
double MomentumMagnitude(double mass, double velocity)
{
  return mass * velocity;
}
double RecoilVelocityMagnitude(double bulletMass, double gunMass, double muzzleVelocity)
{
  return (bulletMass * muzzleVelocity) / gunMass;
}
double ForceMomentumMagnitude(double mass, double velocity1, double velocity2, double time)
{
  return mass * (velocity1 - velocity1) / time;
}
double VelocityMagnitude(double acceleration, double time)
{
  return acceleration * time;
}
double VelocityDisplacementMagnitude(double acceleration, double time)
{
  return (0.5 * acceleration * time * time);
}
double VelocityChange(double acceleration, double displacement)
{
  return Math.sqrt(2 * acceleration * displacement);
}
double ForceMagnitude(double mass, double acceleration)
{
  return mass * acceleration;
}
double SphereVolume(double radius)
{
  return (radius * radius * radius) * 4.0/3.0*Math.PI;
}
double SphereArea(double radius)
{
  return (radius * radius) *Math.PI*4;
}
