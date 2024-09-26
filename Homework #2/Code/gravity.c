#include "gravity.h"

float earth_moon_net_gravity(float earth_r)
{
  const float moon_r = (earth_moon_distance_in_km - earth_r);
  const float earth_gravity = -G * earth_mass_in_kg / (earth_r * earth_r);
  const float moon_gravity = -G * moon_mass_in_kg / (moon_r * moon_r);
  return earth_gravity - moon_gravity;
}
float earth_moon_net_gravity_diff(float earth_r)
{
  const float moon_r = (earth_moon_distance_in_km - earth_r);
  return 2 * G * earth_mass_in_kg / (earth_r * earth_r * earth_r) - 2 * G * moon_mass_in_kg / (moon_r * moon_r * moon_r);
}
void earth_moon_net_gravity_fdf(float x, float *f, float *df)
{
  *f = earth_moon_net_gravity(x);
  *df = earth_moon_net_gravity_diff(x);
}