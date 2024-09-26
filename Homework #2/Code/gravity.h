#ifndef _GRAVITY_H
#define _GRAVITY_H

float earth_moon_net_gravity(float x_from_earth_in_km);
float earth_moon_net_gravity_diff(float x);
void earth_moon_net_gravity_fdf(float x, float *f, float *df);
static const float earth_mass_in_kg = 5.9722e24;
static const float moon_mass_in_kg = 7.34767309e22;
static const float earth_moon_distance_in_km = 384400.0;
static const float G = 6.67430e-11;

#endif /* _GRAVITY_H */