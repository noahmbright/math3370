#include <cmath>

const float g_bar_na = 120;
const float g_bar_k = 36;
const float g_bar_l = 0.3;

const float E_na = 50;
const float E_cl = -40;
const float E_k = -77;
const float E_l = -54.387;

const float c_m = 1.0;
const float phi = 1.0;


inline float alpha_n(float V, float hyperpolarization = 0, float V_half = -55){
	V_half -= hyperpolarization;
	return 0.01 * (V - V_half) / (1 - std::exp( -(V - V_half)/10 ));
}

inline float alpha_m(float V, float hyperpolarization = 0, float V_half = -40){
	V_half -= hyperpolarization;
	return 0.1 * (V - V_half) / (1 - std::exp( -(V - V_half)/10 ));
}

inline float alpha_h(float V, float hyperpolarization = 0, float V_half = -65){
	V_half -= hyperpolarization;
	return 0.07 * std::exp( -(V - V_half)/20 );
}

inline float beta_n(float V, float V_half = -65){
	return 0.125 * std::exp( -(V-V_half)/80);
}

inline float beta_m(float V, float V_half = -65){
	return 4.0 * std::exp(-(V-V_half)/18);
}

inline float beta_h(float V, float V_half = -35){
	return 1.0 / (1 + std::exp( -(V - V_half)/10 ));
}

inline float tau_n(float V, float hyperpolarization = 0){
	return 1.0 / (alpha_n(V, hyperpolarization) + beta_n(V, hyperpolarization));
}

inline float tau_m(float V, float hyperpolarization = 0){
	return 1.0 / (alpha_m(V, hyperpolarization) + beta_m(V, hyperpolarization));
}

inline float tau_h(float V, float hyperpolarization = 0){
	return 1.0 / (alpha_h(V, hyperpolarization) + beta_h(V, hyperpolarization));
}

inline float n_infty(float V, float hyperpolarization = 0){
	return alpha_n(V, hyperpolarization) / (alpha_n(V, hyperpolarization) + beta_n(V, hyperpolarization));
}

inline float m_infty(float V, float hyperpolarization = 0){
	return alpha_m(V, hyperpolarization) / (alpha_m(V, hyperpolarization) + beta_m(V, hyperpolarization));
}

inline float h_infty(float V, float hyperpolarization = 0){
	return alpha_h(V, hyperpolarization) / (alpha_h(V, hyperpolarization) + beta_h(V, hyperpolarization));
}





