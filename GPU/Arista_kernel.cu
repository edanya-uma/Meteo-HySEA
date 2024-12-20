#ifndef _ARISTA_KERNEL_H_
#define _ARISTA_KERNEL_H_

#include <stdio.h>
#include "Matriz.cu"

__device__ double limitador(double p0m, double p0, double p1, double p1p, double S)
{
	double chi;
	double vs;
	double c0,c1;
	double aux;

	c0 = fabs(p1-p0);
	c1 = fabs((p1p-p1)*(S<0.0) + (p0-p0m)*(S>0.0) + 0.5*(p1p-p0m)*(S==0.0));
   
	if (c0 < EPSILON)
		chi = 1.0;
	else {
		aux = fabs((c1-c0)/(c1+c0));
		c1 = 0.5*(c1+c0)*(1.0-aux*aux*aux);
		vs = (c1/(c0+EPSILON))*(c0>c1) + (c0/(c1+EPSILON))*(c1>=c0);
		chi = vs*(1.0+vs)/(1.0+vs*vs);
	}

	return chi;
}

__device__ double limitador_monotono(double p0m, double p0, double p1, double p1p, double S)
{
	double chi;
	double vs;
	double c0,c1;
	double aux,sg;

	c0 = p1-p0;
	c1 = (p1p-p1)*(S<0.0) + (p0-p0m)*(S>0.0) + 0.5*(p1p-p0m)*(S==0.0);
	sg = 1.0*(c0*c1>=0.0);
	c0 = fabs(c0);
	c1 = fabs(c1);
   
	if (c0 < EPSILON)
		chi = 1.0;
	else {
		aux = fabs((c1-c0)/(c1+c0));
		c1 = 0.5*(c1+c0)*(1.0-aux*aux*aux);
		vs = (c1/(c0+EPSILON))*(c0>c1) + (c0/(c1+EPSILON))*(c1>=c0);
		chi = sg*vs*(1.0+vs)/(1.0+vs*vs);
	}

	return chi;
}

__device__ void tratamientoSecoMojado(double h0, double h1, double eta0, double eta1, double H0, double H1,
				double epsilon_h, double *q0n, double *q1n)
{
	if ((h0 < epsilon_h) && (-eta1 > H0))
		*q1n = 0.0;
	if ((h1 < epsilon_h) && (-eta0 > H1))
		*q0n = 0.0;
}

__device__ void positividad(double h0, double h1, double delta_T, double area0, double area1, double *Fmenos)
{
	double dt0;
	double factor = 0.25;
	double alpha;

	if (*Fmenos > 0.0) {
		dt0 = factor*h0*area0/((*Fmenos) + EPSILON);
	}
	else {
		dt0 = factor*h1*area1/(-(*Fmenos) + EPSILON);
	}
	if (delta_T > dt0) {
		alpha = dt0/(delta_T + EPSILON);
		*Fmenos *= alpha;
	}
}

__device__ void positividad_H(double h0, double h1, double cosPhi0, double cosPhi1, double delta_T,
				double area0, double area1, double *Fmenos)
{
	double dt0;
	double factor = 0.25;
	double alpha;

	if (*Fmenos > 0.0) {
		dt0 = factor*h0*cosPhi0*area0/((*Fmenos) + EPSILON);
	}
	else {
		dt0 = factor*h1*cosPhi1*area1/(-(*Fmenos) + EPSILON);
	}
	if (delta_T > dt0) {
		alpha = dt0/(delta_T + EPSILON);
		*Fmenos *= alpha;
	}
}

// ARISTAS VERTICALES

__device__ void procesarAristaVer_2s_1step(double h0, double h1, double q0x, double longitud, double area0,
				double area1, double delta_T, double2 *d_acumulador_1, int pos_vol0, int pos_vol1, double cosPhi0,
				double *Fmas2s, double *Fmenos2s)
{

	*Fmenos2s = q0x*longitud/cosPhi0;

	positividad(h0, h1, delta_T, area0, area1, Fmenos2s);
	*Fmas2s = -(*Fmenos2s);

	d_acumulador_1[pos_vol0].x -= *Fmenos2s;
	d_acumulador_1[pos_vol1].x += *Fmenos2s;
}

__device__ void procesarAristaVer_2s_waf_1step(double h0m, double h0, double h1, double h1p, double q0x, double q1x,
				double H0m, double H0, double H1, double H1p, double normal, double longitud, double area0, double area1,
				double delta_T, double2 *d_acumulador_1, int pos_vol0, int pos_vol1, double epsilon_h, double cosPhi0,
				double hpos, double *Fmas2s, double *Fmenos2s, bool alfa_cero)
{
	double Fmenoswaf;
	double q0n, q1n;
	double eta0m, eta0, eta1, eta1p;
	double hij, uijn;
	double a0, c0;
	double aut1, aut2, max_autovalor;
	double h0rg, h1rg;
	double u0n, u1n;
	double sqrt_h0, sqrt_h1, Hm, deta;
	double alr, blr;
	double sgL, sgR;
	double xiL, xiR;
	double val, alfa, maxheps;
	double ilcos = longitud/cosPhi0;
	double ho0,ho1;
	
	if ((h0 > EPSILON) || (h1 > EPSILON)) {
		q0n = q0x*normal;
		q1n = q1x*normal;
		ho0=h0;
		ho1=h1;
		

		eta0m = h0m - H0m;
		eta0  = h0 - H0;
		eta1  = h1 - H1;
		eta1p = h1p - H1p;
		//tratamientoSecoMojado(h0, h1, eta0, eta1, H0, H1, epsilon_h, &q0n, &q1n);

		maxheps = max(h0,epsilon_h);
		h0rg = M_SQRT2*h0 / sqrt(h0*h0*h0*h0 + maxheps*maxheps*maxheps*maxheps);
		u0n = h0rg*q0n;
		maxheps = max(h1,epsilon_h);
		h1rg = M_SQRT2*h1 / sqrt(h1*h1*h1*h1 + maxheps*maxheps*maxheps*maxheps);
		u1n = h1rg*q1n;
		Hm=min(H0,H1);
		h0=max(h0-H0+Hm,0.0);
		h1=max(h1-H1+Hm,0.0);
		hij = 0.5*(h0 + h1);
		q0n = u0n*h0;
		q1n = u1n*h1;
		
	    *Fmenos2s = ilcos*q0n*normal;

		sqrt_h0 = sqrt(ho0);
		sqrt_h1 = sqrt(ho1);
		uijn = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);

		val = (2.0*ilcos*delta_T) / (area0 + area1);
		c0 = sqrt(0.5*(ho1+ho0));
		aut1 = min(uijn-c0, u0n-sqrt_h0);
		aut2 = max(uijn+c0, u1n+sqrt_h1);
		max_autovalor = max(fabs(aut2), fabs(aut1));

		if (alfa_cero)
			xiL = xiR = 0.0;
		else {
			xiL = limitador_monotono(eta0m,eta0,eta1,eta1p,aut1);
			xiR = limitador_monotono(eta0m,eta0,eta1,eta1p,aut2);
		}

		if (fabs(aut1-aut2) >= EPSILON) {
			sgL = 1.0*(aut1>0.0) - 1.0*(aut1<0.0);
			sgR = 1.0*(aut2>0.0) - 1.0*(aut2<0.0);
			alr = sgL*(1.0-xiL) + val*aut1*xiL;
			blr = sgR*(1.0-xiR) + val*aut2*xiR;
			a0 = (alr-blr)*(aut1*aut2)/(aut2-aut1);
		}
		else {
			a0 = max_autovalor;
		}

		//Hm = min(H0,H1);
		//deta = max(eta1 + Hm, 0.0) - max(eta0 + Hm, 0.0);
		deta=h1-h0;
		Fmenoswaf = 0.5*(q1n + q0n - a0*deta)*ilcos;

		if (alfa_cero)
			alfa = 0.0;
		else
			alfa = (1.0/(1.0 + exp(-15.0*(1.0 - 0.5*hpos/(hij+EPSILON)))))*(hij < hpos) + 1.0*(hij >= hpos);

		*Fmenos2s = alfa*(*Fmenos2s) + (1.0-alfa)*Fmenoswaf;
		positividad(h0, h1, delta_T, area0, area1, Fmenos2s);
		*Fmas2s = -(*Fmenos2s);

		d_acumulador_1[pos_vol0].x -= *Fmenos2s;
		if (pos_vol1 != -1)
			d_acumulador_1[pos_vol1].x += *Fmenos2s;
	}
}

__device__ void procesarAristaVer_2s_2step(TVec3 *W0, TVec3 *W1, double H0, double H1, double2 *d_acumulador_1,
				double2 *d_acumulador_2, int pos_vol0, int pos_vol1, double CFL, double cosPhi0, double longitud,
				double delta_T, double area0, double area1, double cvis, TVec2 *Fmas2s, TVec2 *Fmenos2s, double radio_tierra, double f)
{
	double2 acum;
	double hij, Hm, u0n, u1n;
	double sqrt_h0, sqrt_h1, uijn;
	double max_autovalor, deta;
	double cd, val, DES, a0, c0;
	double ilcos = longitud/cosPhi0;
	double h0,h1,ho0,ho1;
	Hm=min(H0,H1);

	ho0=h0=W0->x;
	ho1=h1=W1->x;
	u0n = W0->y/W0->x;
	u1n = W1->y/W1->x;
	h0=max(h0-H0+Hm,0.0);
	h1=max(h1-H1+Hm,0.0);
	sqrt_h0 = sqrt(ho0);
	sqrt_h1 = sqrt(ho1);
	hij=0.5*(h0+h1);
	uijn = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);
	c0 = sqrt(0.5*(ho0+ho1));
	max_autovalor = fabs(uijn) + c0;
	deta=h1-h0;
	//deta = max(W1->x - H1 + Hm, 0.0) - max(W0->x - H0 + Hm, 0.0);

	//NUEVO
	double dx;
	dx = area0/(longitud);
	//-----------
	double u0t=W0->z/W0->x;
	double u1t=W1->z/W1->x;
	double q0n=u0n*h0;
	double q1n=u1n*h1;
	if (uijn < 0.0) {
		Fmenos2s->x = u1n*(q1n)*ilcos;
		Fmas2s->x = -Fmenos2s->x;
		//Fmenos2s->x += hij*deta*ilcos;
		//Fmenos2s->x += hij*(deta-f*cosPhi0*0.5*(u0t+u1t)*(dx))*ilcos; //NUEVO
		Fmenos2s->x += hij*(deta-f*cosPhi0*(u0t)*(dx))*ilcos; //NUEVO
		Fmenos2s->y = u1n*(u1t*h1)*ilcos;
		Fmas2s->y = -Fmenos2s->y;
	}
	else {
		Fmenos2s->x = u0n*(q0n)*ilcos;
		Fmas2s->x = -Fmenos2s->x;
		//Fmenos2s->x += hij*deta*ilcos;
		//Fmenos2s->x += hij*(deta-radio_tierra*f*cosPhi0*cosPhi0*0.5*(u0t+u1t)*(dx))*ilcos;	//NUEVO
		//Fmenos2s->x += hij*(deta-f*cosPhi0*0.5*(u0t+u1t)*(dx))*ilcos; //NUEVO
		Fmenos2s->x += hij*(deta-f*cosPhi0*(u0t)*(dx))*ilcos; //NUEVO
		Fmenos2s->y = u0n*(u0t*h0)*ilcos;
		Fmas2s->y = -Fmenos2s->y;
	}

	cd = 0.5*(1.0-cvis);
	//cd=0.5;

	val = (2.0*ilcos*delta_T) / (area0 + area1);
	DES = cd*val*max_autovalor*max_autovalor*(q1n-q0n)*ilcos;
	Fmenos2s->x -= DES;
	Fmas2s->x += DES;

	a0 = ilcos*max_autovalor;
	d_acumulador_1[pos_vol0].y += a0;
	acum = d_acumulador_2[pos_vol0];
	acum.x -= Fmenos2s->x;
	acum.y -= Fmenos2s->y;
	d_acumulador_2[pos_vol0] = acum;

	d_acumulador_1[pos_vol1].y += a0;
	acum = d_acumulador_2[pos_vol1];
	acum.x -= Fmas2s->x;
	acum.y -= Fmas2s->y;
	d_acumulador_2[pos_vol1] = acum;
}

__device__ void procesarAristaVer_2s_waf_2step(TVec3 *W0m, TVec3 *W0, TVec3 *W1, TVec3 *W1p, double H0m,
				double H0, double H1, double H1p, double normal, double longitud, double area0, double area1,
				double delta_T, double2 *d_acumulador_1, double2 *d_acumulador_2, int pos_vol0, int pos_vol1,
				double CFL, double epsilon_h, double cosPhi0, double cvis, double hpos, TVec2 *Fmas2,
				TVec2 *Fmenos2, bool alfa_cero, double radio_tierra, double f)
{
	TVec2 DES;
	TVec3 Fmas3, Fmenos3;
	TVec3 Fmas2s, Fmenos2s;
	TVec2 fl0, fl1;
	double2 acum;
	double hij, uijn, uijt,h0,h1;
	double a0, a1, c0;
	double aut1, aut2, max_autovalor;
	double h0rg, h1rg;
	double u0n, u1n;
	double q0mt, q0t, q1t, q1pt;
	double sqrt_h0, sqrt_h1, Hm, deta, dq;
	TVec2 W0_rot, W1_rot;
	double eta0m, eta0, eta1, eta1p;
	double h0m, h1p;
	double p0m, p0, p1, p1p;
	double alr, blr;
	double sgL, sgR;
	double xiL, xiR;
	double val, aux, maxheps;
	double cd, alfa, des;
	double ilcos = longitud/cosPhi0;
	double ho0,ho1;

	if ((W0->x > EPSILON) || (W1->x > EPSILON)) {
		ho0=h0 = W0_rot.x = W0->x;
		W0_rot.y = W0->y*normal;
		ho1=h1 = W1_rot.x = W1->x;
		W1_rot.y = W1->y*normal;
		//hij = 0.5*(h0 + h1);
		h1p = W1p->x;
		h0m = W0m->x;
		eta0m = h0m - H0m;
		eta0  = h0 - H0;
		eta1  = h1 - H1;
		eta1p = h1p - H1p;
		Hm=min(H0,H1);

		//tratamientoSecoMojado(h0, h1, eta0, eta1, H0, H1, epsilon_h, &(W0_rot.y), &(W1_rot.y)); //pensar si modificar quitandole a H0,H1 la presion

		maxheps = max(h0,epsilon_h);
		h0rg = M_SQRT2*h0 / sqrt(h0*h0*h0*h0 + maxheps*maxheps*maxheps*maxheps);
		u0n = h0rg*W0_rot.y;
		maxheps = max(h1,epsilon_h);
		h1rg = M_SQRT2*h1 / sqrt(h1*h1*h1*h1 + maxheps*maxheps*maxheps*maxheps);
		u1n = h1rg*W1_rot.y;
		W0_rot.x=h0=max(h0-H0+Hm,0.0);
		W1_rot.x=h1=max(h1-H1+Hm,0.0);
		hij=0.5*(h0+h1);		

		
		W0_rot.y = u0n*h0;
		W1_rot.y = u1n*h1;

		sqrt_h0 = sqrt(ho0);
		sqrt_h1 = sqrt(ho1);
		uijn = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);

		val = (2.0*ilcos*delta_T) / (area0 + area1);
		c0 = sqrt(0.5*(ho0+ho1));
		aut1 = min(uijn-c0, u0n-sqrt_h0);
		aut2 = max(uijn+c0, u1n+sqrt_h1);
		max_autovalor = max(fabs(aut2), fabs(aut1));

		fl1.x = W1_rot.y;
		fl1.y = u1n*W1_rot.y;
		fl0.x = W0_rot.y;
		fl0.y = u0n*W0_rot.y;
		dq = fl1.x-fl0.x;

		//Hm = min(H0,H1);
		//deta = max(eta1 + Hm, 0.0) - max(eta0 + Hm, 0.0);
		deta=h1-h0;
		//NUEVO
		double v0 = W0->z*h0rg;
		double v1 = W1->z*h1rg; 
		double dx;
		//dx = area0/(radio_tierra*longitud);
		dx = area0/(longitud);
		//-----------

		

		if (uijn < 0.0) {
			Fmenos2s.y = u1n*W1_rot.y*ilcos*normal;
			Fmas2s.y = -Fmenos2s.y;
			//Fmenos2s.y += hij*deta*ilcos;
			//Fmenos2s.y += hij*(deta-radio_tierra*f*cosPhi0*cosPhi0*0.5*(v0+v1)*dx)*ilcos;	//NUEVO
			//Fmenos2s.y += hij*(deta-f*cosPhi0*0.5*(v0+v1)*dx)*ilcos;	//NUEVO
			Fmenos2s.y += hij*(deta-f*cosPhi0*(v0)*dx)*ilcos;	//NUEVO
			Fmenos2s.z = u1n*v1*h1*ilcos;
			Fmas2s.z = -Fmenos2s.z;
		}
		else {
			Fmenos2s.y = u0n*W0_rot.y*ilcos*normal;
			Fmas2s.y = -Fmenos2s.y;
			//Fmenos2s.y += hij*deta*ilcos;
			//Fmenos2s.y += hij*(deta-radio_tierra*f*cosPhi0*cosPhi0*0.5*(v0+v1)*dx)*ilcos;	//NUEVO
			//Fmenos2s.y += hij*(deta-f*cosPhi0*0.5*(v0+v1)*dx)*ilcos;	//NUEVO
			Fmenos2s.y += hij*(deta-f*cosPhi0*(v0)*dx)*ilcos;	//NUEVO
			Fmenos2s.z = u0n*v0*h0*ilcos;
			Fmas2s.z = -Fmenos2s.z;
		}

		cd = 0.5*(1.0 - cvis/(1.0 + exp(-20.0*(1.0 - hpos/hij))));
		//cd=0.5;
		des = cd*val*max_autovalor*max_autovalor*dq*ilcos*normal;
		Fmenos2s.y -= des;
		Fmas2s.y += des;

		if (alfa_cero)
			xiL = xiR = 0.0;
		else {
			xiL = limitador_monotono(eta0m,eta0,eta1,eta1p,aut1);
			xiR = limitador_monotono(eta0m,eta0,eta1,eta1p,aut2);
		}

		if (fabs(aut1-aut2) >= EPSILON) {
			sgL = 1.0*(aut1>0.0) - 1.0*(aut1<0.0);
			sgR = 1.0*(aut2>0.0) - 1.0*(aut2<0.0);
			alr = sgL*(1.0-xiL) + val*aut1*xiL;
			blr = sgR*(1.0-xiR) + val*aut2*xiR;
			aux = 1.0/(aut2-aut1);
			a0 = (alr-blr)*(aut1*aut2)*aux;
			a1 = (blr*aut2-alr*aut1)*aux;
		}
		else {
			a0 = max_autovalor;
			a1 = 0.0;
		}

		Fmenos2->x = Fmas2->x = 0.5*dq;
		//Fmenos2->y = Fmas2->y = 0.5*(fl1.y-fl0.y + hij*deta);
		Fmenos2->y = Fmas2->y = 0.5*(fl1.y-fl0.y + hij*(deta-f*cosPhi0*0.5*(v0+v1)*dx));	//NUEVO
		//DES.x = 0.5*a0*deta;
		DES.x = 0.5*a0*(deta-f*cosPhi0*0.5*(v0+v1)*dx)+a1*Fmenos2->x; //NUEVO
		DES.y = 0.5*a0*dq + a1*Fmenos2->y;
		Fmenos2->x += fl0.x-DES.x;
		Fmenos2->y += fl0.y-DES.y;
		Fmas2->x += DES.x-fl1.x;
		Fmas2->y += DES.y-fl1.y;

		Fmenos2->x *= ilcos;
		Fmenos2->y *= ilcos;
		Fmas2->x *= ilcos;
		Fmas2->y *= ilcos;

		q0mt = W0m->z*normal;
		q0t = W0->z*normal;
		q1t = W1->z*normal;
		q1pt = W1p->z*normal;
		maxheps = max(h0m,epsilon_h);
		p0m = M_SQRT2*h0m*q0mt / sqrt(h0m*h0m*h0m*h0m + maxheps*maxheps*maxheps*maxheps);
		p0 = u0n = h0rg*q0t;
		p1 = u1n = h1rg*q1t;
		maxheps = max(h1p,epsilon_h);
		p1p = M_SQRT2*h1p*q1pt / sqrt(h1p*h1p*h1p*h1p + maxheps*maxheps*maxheps*maxheps);

		a0 = -1.0*(Fmenos2->x<0.0) + 1.0*(Fmenos2->x>0.0);
		uijn = fabs(uijn)*a0;
		if (alfa_cero)
			xiL = 0.0;
		else
			xiL = limitador_monotono(p0m, p0, p1, p1p, uijn);
		uijt = 0.5*(u0n+u1n) - 0.5*((1.0-xiL)*a0 + xiL*uijn*val)*(u1n-u0n);

		aux = Fmenos2->x*uijt;
		Fmas3.y = Fmas2->y*normal;
		Fmas3.z = -aux*normal;
		Fmenos3.y = Fmenos2->y*normal;
		Fmenos3.z = aux*normal;

		if (alfa_cero)
			alfa = 0.0;
		else
			alfa = (1.0/(1.0 + exp(-15.0*(1.0 - 0.5*hpos/(hij+EPSILON)))))*(hij < hpos) + 1.0*(hij >= hpos);

		Fmenos2->x = alfa*Fmenos2s.y + (1.0-alfa)*Fmenos3.y;
		Fmenos2->y = alfa*Fmenos2s.z + (1.0-alfa)*Fmenos3.z;

		Fmas2->x = alfa*Fmas2s.y + (1.0-alfa)*Fmas3.y;
		Fmas2->y = alfa*Fmas2s.z + (1.0-alfa)*Fmas3.z;

		if (max_autovalor < epsilon_h)
			max_autovalor += epsilon_h;

		a0 = ilcos*max_autovalor;
		d_acumulador_1[pos_vol0].y += a0;
		acum = d_acumulador_2[pos_vol0];
		acum.x -= Fmenos2->x;
		acum.y -= Fmenos2->y;
		d_acumulador_2[pos_vol0] = acum;

		if (pos_vol1 != -1) {
			d_acumulador_1[pos_vol1].y += a0;
			acum = d_acumulador_2[pos_vol1];
			acum.x -= Fmas2->x;
			acum.y -= Fmas2->y;
			d_acumulador_2[pos_vol1] = acum;
		}
	}
}

__global__ void procesarAristasVerNivel0Paso1GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double *d_altoVolumenes, int num_volx, int num_voly,
				double borde1, double borde2, double delta_T, double2 *d_acumulador_1,
				double CFL, double epsilon_h, double hpos, double cvis, int tipo)
{
	double2 datos_vol;
	double normal;
	double Fmas, Fmenos;
	double area0, area1, longitud;
	int pos_x_hebra, pos_y_hebra;
	int pos, pos_vol0, pos_vol1;
	int ppv, ppv0, front;
	double H0m, H0, H1, H1p, aux;
	double h0m, h0, h1, h1p;
	double q0x, q1x;
	double cosPhi0;
	bool alfa_cero;

	pos_x_hebra = 2*(blockIdx.x*NUM_HEBRASX_ARI + threadIdx.x);
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_ARI + threadIdx.y;
	if (tipo == 2) pos_x_hebra++;

	if ((pos_x_hebra <= num_volx) && (pos_y_hebra < num_voly)) {
		pos_vol0 = pos_y_hebra*num_volx + pos_x_hebra-1;
		pos_vol1 = pos_vol0 + 1;
		longitud = d_altoVolumenes[pos_y_hebra];
		datos_vol = d_datosVolumenes_3[pos_y_hebra];
		area0 = datos_vol.x;
		cosPhi0 = datos_vol.y;
		area1 = area0;
		normal = 1.0;
		front = 0;

		//Primer ppv (h0)
		ppv = pos_x_hebra-1;
		ppv0 = ppv*(ppv >= 0);
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		h0 = datos_vol.x;
		H0 = datos_vol.y;
		q0x = d_datosVolumenes_2[pos].x;
		if (ppv < 0) {
			//aqui creo que habria que hacer algo si la frontera es periodica
			if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
				pos_vol0 = pos_y_hebra*num_volx;
				//pos_vol1 = pos_y_hebra*num_volx + num_volx-1;
				pos_vol1 = -1.0;
				normal = -1.0;
			}
			else{
				q0x *= borde1;
				pos_vol0 = pos_y_hebra*num_volx;
				pos_vol1 = -1;
				normal = -1.0;
				front = 1;
			}
			
		}

		//Segundo ppv (h1)
		ppv = pos_x_hebra;
		if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
			ppv0 = (num_volx-1)*(ppv==0)+0*(ppv==num_volx)+ppv*(0<ppv && ppv<num_volx);
		}
		else{
			ppv0 = ppv*(ppv < num_volx) + (num_volx-1)*(ppv >= num_volx);
		}
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		h1 = datos_vol.x;
		H1 = datos_vol.y;
		q1x = d_datosVolumenes_2[pos].x;
		if (ppv >= num_volx) {
			if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
				pos_vol0 = pos_y_hebra*num_volx + num_volx-1;
				pos_vol1 = -1;
			}
			else{
				q1x *= borde2;
				pos_vol0 = (pos_y_hebra+1)*num_volx - 1;
				pos_vol1 = -1; 
			}
			
		}

		//Tercer ppv (h1m)
		ppv = pos_x_hebra-2;
		if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
			ppv0 = (num_volx-2)*(ppv == -2)+(num_volx-1)*(ppv == -1)+ppv*(ppv >= 0);
		}
		else{
			ppv0 = ppv*(ppv >= 0);
		}
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		h0m = datos_vol.x;
		H0m = datos_vol.y;

		//Cuarto ppv (h1p)
		ppv = pos_x_hebra+1;
		if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
			ppv0 = 0*(ppv == num_volx)+1*(ppv == (num_volx+1))+ppv*(ppv < num_volx);
		}
		else{
			ppv0 = ppv*(ppv < num_volx) + (num_volx-1)*(ppv >= num_volx);
		}
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		h1p = datos_vol.x;
		H1p = datos_vol.y;

		alfa_cero = ( ((pos_x_hebra < 2) || (pos_x_hebra > num_volx-2)) ? true : false);
		//alfa_cero = ( ((pos_x_hebra < 2) || (pos_x_hebra > num_volx-2)|| (cosPhi0<0.09)) ? true : false);
		//alfa_cero = true;
		if ((h0 > hpos) && (h1 > hpos) && (h0m > hpos) && (h1p > hpos) && (! alfa_cero)) {
            procesarAristaVer_2s_1step(h0, h1, q0x, longitud, area0, area1, delta_T,
				d_acumulador_1, pos_vol0, pos_vol1, cosPhi0, &Fmas, &Fmenos);
        }
        else {
            if (front == 1) {
                aux = H0m; H0m = H1p; H1p = aux;
                aux = H0;  H0 = H1;   H1 = aux;
                aux = h0m; h0m = h1p; h1p = aux;
                aux = h0;  h0 = h1;   h1 = aux;
                aux = q0x; q0x = q1x; q1x = aux;
            }
            procesarAristaVer_2s_waf_1step(h0m, h0, h1, h1p, q0x, q1x, H0m, H0, H1, H1p, normal,
				longitud, area0, area1, delta_T, d_acumulador_1, pos_vol0, pos_vol1, epsilon_h,
				cosPhi0, hpos, &Fmas, &Fmenos, alfa_cero);
        }

		
    }
}

__global__ void procesarAristasVerNivel0Paso2GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double *d_altoVolumenes, int num_volx, int num_voly, double borde1,
				double borde2, double delta_T, double2 *d_acumulador_1, double2 *d_acumulador_2,
				double CFL, double epsilon_h, double hpos, double cvis, int tipo, double radio_tierra, double *sinPhi, double T)	//NUEVO
{
	double2 datos_vol;
	double normal;
	TVec2 Fmas2, Fmenos2;
	double area0, area1, longitud;
	int pos_x_hebra, pos_y_hebra;
	int pos, pos_vol0, pos_vol1;
	int ppv, ppv0, front;
	TVec3 W0, W1, W0m, W1p, vaux;
	double H0m, H0, H1, H1p, aux;
	double cosPhi0;
	bool alfa_cero;
	double f; //NUEVO

	pos_x_hebra = 2*(blockIdx.x*NUM_HEBRASX_ARI + threadIdx.x);
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_ARI + threadIdx.y;
	if (tipo == 2) pos_x_hebra++;

	if ((pos_x_hebra <= num_volx) && (pos_y_hebra < num_voly)) {
		pos_vol0 = pos_y_hebra*num_volx + pos_x_hebra-1;
		pos_vol1 = pos_vol0 + 1;
		longitud = d_altoVolumenes[pos_y_hebra];
		datos_vol = d_datosVolumenes_3[pos_y_hebra];
		area0 = datos_vol.x;
		cosPhi0 = datos_vol.y;
		area1 = area0;
		front = 0;
		normal = 1.0;
		f = 2.0*EARTH_ANGULAR_VELOCITY*sinPhi[pos_y_hebra]*T;	//NUEVO (multiplico por T para normalizar el parametro de coriolis)

		//Primer ppv (W0)
		ppv = pos_x_hebra-1;
		ppv0 = ppv*(ppv >= 0);
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		W0.x = datos_vol.x;
		H0 = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W0.y = datos_vol.x;
		W0.z = datos_vol.y;
		if (ppv < 0) {
			if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
				pos_vol0 = pos_y_hebra*num_volx;
				pos_vol1 = -1;
				normal = -1.0;
			}
			else{
				W0.y *= borde1;
				pos_vol0 = pos_y_hebra*num_volx;
				pos_vol1 = -1;
				front = 1;
				normal = -1.0;
			}
			
		}

		//Segundo ppv (W1)
		ppv = pos_x_hebra;
		if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
			ppv0 = (num_volx-1)*(ppv==0)+0*(ppv==num_volx)+ppv*(0<ppv && ppv<num_volx);
		}
		else{
			ppv0 = ppv*(ppv < num_volx) + (num_volx-1)*(ppv >= num_volx);
		}
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		W1.x = datos_vol.x;
		H1 = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W1.y = datos_vol.x;
		W1.z = datos_vol.y;
		if (ppv >= num_volx) {
			if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
				pos_vol0 = pos_y_hebra*num_volx + num_volx-1;
				pos_vol1 = -1;
			}
			else{
				W1.y *= borde2;
				pos_vol0 = (pos_y_hebra+1)*num_volx - 1;
				pos_vol1 = -1; 
			}
		}

		//Tercer ppv (W0m)
		ppv = pos_x_hebra-2;
		if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
			ppv0 = (num_volx-2)*(ppv == -2)+(num_volx-1)*(ppv == -1)+ppv*(ppv >= 0);
		}
		else{
			ppv0 = ppv*(ppv >= 0);
		}
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		W0m.x = datos_vol.x;
		H0m = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W0m.y = datos_vol.x;
		W0m.z = datos_vol.y;

		//Cuarto ppv (W1p)
		ppv = pos_x_hebra+1;
		if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
			ppv0 = 0*(ppv == num_volx)+1*(ppv == (num_volx+1))+ppv*(ppv < num_volx);
		}
		else{
			ppv0 = ppv*(ppv < num_volx) + (num_volx-1)*(ppv >= num_volx);
		}
		pos = pos_y_hebra*num_volx + ppv0;
		datos_vol = d_datosVolumenes_1[pos];
		W1p.x = datos_vol.x;
		H1p = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W1p.y = datos_vol.x;
		W1p.z = datos_vol.y;

		alfa_cero = ( ((pos_x_hebra < 2) || (pos_x_hebra > num_volx-2)) ? true : false);
		//alfa_cero = ( ((pos_x_hebra < 2) || (pos_x_hebra > num_volx-2)|| (cosPhi0<0.09)) ? true : false);
		//alfa_cero = true;
		
		if ((W0.x > hpos) && (W1.x > hpos) && (W0m.x > hpos) && (W1p.x > hpos) && (! alfa_cero)) {
            procesarAristaVer_2s_2step(&W0, &W1, H0, H1, d_acumulador_1, d_acumulador_2, pos_vol0,
				pos_vol1, CFL, cosPhi0, longitud, delta_T, area0, area1, cvis, &Fmas2, &Fmenos2, radio_tierra, f);	//NUEVO (MODIFICADO)
        }
        else {
            if (front == 1) {
                v_copy3(&W0, &vaux);  v_copy3(&W1, &W0);  v_copy3(&vaux, &W1);
                v_copy3(&W0m, &vaux);  v_copy3(&W1p, &W0m);  v_copy3(&vaux, &W1p);
                aux = H0m; H0m = H1p; H1p = aux;
                aux = H0;  H0 = H1;   H1 = aux;
            }

            procesarAristaVer_2s_waf_2step(&W0m, &W0, &W1, &W1p, H0m, H0, H1, H1p, normal, longitud,
				area0, area1, delta_T, d_acumulador_1, d_acumulador_2, pos_vol0, pos_vol1, CFL,
				epsilon_h, cosPhi0, cvis, hpos, &Fmas2, &Fmenos2, alfa_cero, radio_tierra, f);	//NUEVO (MODIFICADO)
        }

		
	}
}

// ARISTAS HORIZONTALES

__device__ void procesarAristaHor_2s_1step(double h0, double h1, double q0y, double longitud, double area0,
				double area1, double delta_T, double2 *d_acumulador_1, int pos_vol0, int pos_vol1, double cosPhi0,
				double cosPhi1, double *Fmas2s, double *Fmenos2s)
{
	*Fmenos2s = longitud*q0y*cosPhi0;
	positividad_H(h0, h1, cosPhi0, cosPhi1, delta_T, area0, area1, Fmenos2s);
	*Fmas2s = -(*Fmenos2s);

	if (pos_vol0 >= 0)
		d_acumulador_1[pos_vol0].x -= (*Fmenos2s)/cosPhi0;
	d_acumulador_1[pos_vol1].x += (*Fmenos2s)/cosPhi1;
}

__device__ void procesarAristaHor_2s_waf_1step(double h0m, double h0, double h1, double h1p, double q0y, double q1y,
				double H0m, double H0, double H1, double H1p, double normal, double longitud, double area0, double area1,
				double delta_T, double2 *d_acumulador_1, int pos_vol0, int pos_vol1, double epsilon_h, double cosPhi0,
				double cosPhi1, double hpos, double *Fmas2s, double *Fmenos2s, bool alfa_cero)
{
	double Fmenoswaf;
	double hij, uijn;
	double fl0, fl1;
	double q0n, q1n;
	double a0, c0;
	double aut1, aut2, max_autovalor;
	double h0rg, h1rg;
	double u0n, u1n;
	double sqrt_h0, sqrt_h1, deta, Hm;
	double eta0m, eta0, eta1, eta1p;
	double alr, blr;
	double sgL, sgR;
	double xiL, xiR;
	double val, alfa, maxheps;
	double ho0,ho1;

	if ((h0 > EPSILON) || (h1 > EPSILON)) {
		q0n = q0y*normal;
		q1n = q1y*normal;
		ho0=h0;
		ho1=h1;
		eta0m = h0m - H0m;
		eta0  = h0 - H0;
		eta1  = h1 - H1;
		eta1p = h1p - H1p;

		//tratamientoSecoMojado(h0, h1, eta0, eta1, H0, H1, epsilon_h, &q0n, &q1n);

		maxheps = max(h0,epsilon_h);
		h0rg = M_SQRT2*h0 / sqrt(h0*h0*h0*h0 + maxheps*maxheps*maxheps*maxheps);
		u0n = h0rg*q0n;
		maxheps = max(h1,epsilon_h);
		h1rg = M_SQRT2*h1 / sqrt(h1*h1*h1*h1 + maxheps*maxheps*maxheps*maxheps);
		u1n = h1rg*q1n;
		Hm=min(H0,H1);
		h0=max(h0-H0+Hm,0.0);
		h1=max(h1-H1+Hm,0.0);
		q0n = u0n*h0;
		q1n = u1n*h1;
		hij = 0.5*(h0 + h1);
		c0 = sqrt(0.5*(ho0+ho1));
		*Fmenos2s = longitud*q0n*cosPhi0*normal;

		sqrt_h0 = sqrt(ho0);
		sqrt_h1 = sqrt(ho1);
		uijn = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);

		val = (2.0*longitud*delta_T) / (area0 + area1);
		aut1 = min(uijn-c0, u0n - sqrt_h0);
		aut2 = max(uijn+c0, u1n + sqrt_h1);
		max_autovalor = max(fabs(aut2), fabs(aut1));

		if (alfa_cero)
			xiL = xiR = 0.0;
		else {
			xiL = limitador_monotono(eta0m,eta0,eta1,eta1p,aut1);
			xiR = limitador_monotono(eta0m,eta0,eta1,eta1p,aut2);
		}

		if (fabs(aut1-aut2) >= EPSILON) {
			sgL = 1.0*(aut1>0.0) - 1.0*(aut1<0.0);
			sgR = 1.0*(aut2>0.0) - 1.0*(aut2<0.0);
			alr = sgL*(1.0-xiL) + val*aut1*xiL;
			blr = sgR*(1.0-xiR) + val*aut2*xiR;
			a0 = (alr-blr)*(aut1*aut2)/(aut2-aut1);
		}
		else {
			a0 = max_autovalor;
		}
		fl1 = q1n*cosPhi1;
		fl0 = q0n*cosPhi0;

		//Hm = min(H0,H1);
		//deta = 0.5*(cosPhi0+cosPhi1)*(max(eta1+Hm, 0.0) - max(eta0+Hm, 0.0)) - uijn*uijn*(cosPhi1-cosPhi0);
		deta=0.5*(cosPhi0+cosPhi1)*(h1-h0);
		Fmenoswaf = 0.5*(fl0+fl1-a0*deta)*longitud;

		if (alfa_cero)
			alfa = 0.0;
		else
			alfa = (1.0/(1.0 + exp(-15.0*(1.0 - 0.5*hpos/(hij+EPSILON)))))*(hij < hpos) + 1.0*(hij >= hpos);

		*Fmenos2s = alfa*(*Fmenos2s) + (1.0-alfa)*Fmenoswaf;
		positividad_H(h0, h1, cosPhi0, cosPhi1, delta_T, area0, area1, Fmenos2s);
		*Fmas2s = -(*Fmenos2s);

		if (pos_vol0 >= 0)
			d_acumulador_1[pos_vol0].x -= (*Fmenos2s)/cosPhi0;
		if (pos_vol1 != -1)
			d_acumulador_1[pos_vol1].x += (*Fmenos2s)/cosPhi1;
	}
}

__device__ void procesarAristaHor_2s_2step(TVec3 *W0, TVec3 *W1, double H0, double H1, double longitud,
				double area0, double area1, double delta_T, double2 *d_acumulador_1, double2 *d_acumulador_2,
				int pos_vol0, int pos_vol1, double CFL, double epsilon_h, double cosPhi0, double cosPhi1,
				double cvis, double hpos, TVec2 *Fmas2s, TVec2 *Fmenos2s, double radio_tierra, double f)
{
	double2 acum;
	double a0;
	double hij, Hm, u0n, u1n,u0t,u1t;
	double sqrt_h0, sqrt_h1, uijn;
	double max_autovalor, deta, c0;
	double val, cd, DES;
	double icos0 = 1.0/cosPhi0;
	double icos1 = 1.0/cosPhi1;
	double h0,h1,q1n,q0n;
	double ho0,ho1;
	Hm=min(H0,H1);
	ho0=W0->x;
	ho1=W1->x;
	h0=max(W0->x-H0+Hm,0.0);
	h1=max(W1->x-H1+Hm,0.0);
	hij = 0.5*(h0+h1);
	u0n = W0->z/W0->x;
	u1n = W1->z/W1->x;
	u0t = W0->y/W0->x;
	u1t = W1->y/W1->x;
	q0n=u0n*h0;
	q1n=u1n*h1;

	sqrt_h0 = sqrt(ho0);
	sqrt_h1 = sqrt(ho1);
	uijn = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);
	c0 = sqrt(0.5*(ho0+ho1));
	max_autovalor = fabs(uijn) + c0;
	//deta = 0.5*(cosPhi0+cosPhi1)*(h1-h0)- uijn*uijn*(cosPhi1-cosPhi0);
	deta = 0.5*(cosPhi0+cosPhi1)*(h1-h0);
	//NUEVO
	double dy;
	dy = area0/(longitud);
	//------------------

	if (uijn < 0.0) {
		Fmenos2s->y = u1n*(q1n)*longitud;
		Fmas2s->y = -Fmenos2s->y;
		Fmenos2s->y *= cosPhi1*icos0;
		//Fmenos2s->y += hij*deta*longitud*icos0;
		//Fmenos2s->y += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*0.5*(u0t+u1t)*(dy))*longitud*icos0;	//NUEVO
		Fmenos2s->y += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*(u0t)*(dy))*longitud*icos0;	//NUEVO
		Fmenos2s->x = u1n*(u1t*h1)*longitud;
		Fmas2s->x = -Fmenos2s->x;
	}
	else {
		Fmenos2s->y = u0n*(q0n)*longitud;
		Fmas2s->y = -Fmenos2s->y*cosPhi0*icos1;
		//Fmenos2s->y += hij*deta*longitud*icos0;
		//Fmenos2s->y += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*0.5*(u0t+u1t)*(dy))*longitud*icos0;	//NUEVO
		Fmenos2s->y += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*(u0t)*(dy))*longitud*icos0;	//NUEVO
		Fmenos2s->x = u0n*(u0t*h0)*longitud;
		Fmas2s->x = -Fmenos2s->x;
	}

	val = (2.0*longitud*delta_T) / (area0 + area1);
	cd = 0.5*(1.0-cvis);
	//cd=0.5;
	DES = cd*val*max_autovalor*max_autovalor*(q1n*cosPhi1-q0n*cosPhi0)*longitud;
	Fmenos2s->y -= DES*icos0;
	Fmas2s->y += DES*icos1;

	a0 = longitud*max_autovalor;
	if (pos_vol0 >= 0) {
		d_acumulador_1[pos_vol0].y += a0;
		acum = d_acumulador_2[pos_vol0];
		acum.x -= Fmenos2s->x;
		acum.y -= Fmenos2s->y;
		d_acumulador_2[pos_vol0] = acum;
	}

	d_acumulador_1[pos_vol1].y += a0;
	acum = d_acumulador_2[pos_vol1];
	acum.x -= Fmas2s->x;
	acum.y -= Fmas2s->y;
	d_acumulador_2[pos_vol1] = acum;
}

__device__ void procesarAristaHor_2s_waf_2step(TVec3 *W0m, TVec3 *W0, TVec3 *W1, TVec3 *W1p, double H0m,
				double H0, double H1, double H1p, double normal, double longitud, double area0, double area1,
				double delta_T, double2 *d_acumulador_1, double2 *d_acumulador_2, int pos_vol0, int pos_vol1,
				double CFL, double epsilon_h, double cosPhi0, double cosPhi1, double cvis, double hpos,
				TVec2 *Fmas2, TVec2 *Fmenos2, bool alfa_cero, double radio_tierra, double f)
{
	TVec2 DES;
	TVec3 Fmas3, Fmenos3;
	TVec3 Fmas2s, Fmenos2s;
	TVec2 fl0, fl1;
	double2 acum;
	double hij, uijn, uijt;
	double a0, a1, c0;
	double aut1, aut2, max_autovalor;
	double h0rg, h1rg;
	double u0n, u1n;
	double q0mt, q0t, q1t, q1pt;
	double h0m, h0, h1, h1p, ho0, ho1;
	double sqrt_h0, sqrt_h1, deta, dq, Hm;
	TVec2 W0_rot, W1_rot;
	double p0m, p0, p1, p1p;
	double alr, blr;
	double sgL, sgR;
	double xiL, xiR;
	double val, aux, maxheps;
	double cd, alfa, des;
	double icos0,icos1;

    if ((W0->x > EPSILON) || (W1->x > EPSILON)) {
        icos0 = 1.0/cosPhi0;
        icos1 = 1.0/cosPhi1;
        ho0=h0=W0_rot.x = W0->x;
        W0_rot.y = W0->z*normal;
        ho1=h1=W1_rot.x = W1->x;
        W1_rot.y = W1->z*normal;
        

        p0m = W0m->x - H0m;
        p0  = W0->x - H0;
        p1  = W1->x - H1;
        p1p = W1p->x - H1p;
        //tratamientoSecoMojado(h0, h1, p0, p1, H0, H1, epsilon_h, &(W0_rot.y), &(W1_rot.y)); //pensar si modificar

		maxheps = max(h0,epsilon_h);
        h0rg = M_SQRT2*h0 / sqrt(h0*h0*h0*h0 + maxheps*maxheps*maxheps*maxheps);
        u0n = h0rg*W0_rot.y;
		maxheps = max(h1,epsilon_h);
        h1rg = M_SQRT2*h1 / sqrt(h1*h1*h1*h1 + maxheps*maxheps*maxheps*maxheps);
        u1n = h1rg*W1_rot.y;
		Hm = min(H0,H1);
		h0=max(h0-H0+Hm,0.0);
		h1=max(h1-H1+Hm,0.0);

		W0_rot.x=h0;
        W1_rot.x=h1;
        hij = 0.5*(h0 + h1);
        c0 = sqrt(0.5*(ho0+ho1));
        sqrt_h0 = sqrt(ho0);
        sqrt_h1 = sqrt(ho1);
        uijn = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);
        val = (2.0*longitud*delta_T) / (area0 + area1);
        aut1 = min(uijn-c0, u0n - sqrt_h0);
        aut2 = max(uijn+c0, u1n + sqrt_h1);
        max_autovalor = max(fabs(aut2), fabs(aut1));
		double deta0=0.5*(cosPhi0+cosPhi1)*(h1-h0);
        deta=deta0 - uijn*uijn*(cosPhi1-cosPhi0);
        //deta = 0.5*(cosPhi0+cosPhi1)*(max(p1+Hm, 0.0) - max(p0+Hm, 0.0)) - uijn*uijn*(cosPhi1-cosPhi0);

		//NUEVO
		double u0, u1;
		u0 = W0->y*h0rg;
		u1 = W1->y*h1rg;

		double dy;
		dy = area0/(longitud);
		//------------------
        
        if (uijn < 0.0) {
            Fmenos2s.z = u1n*W1_rot.y*longitud*normal;
            Fmas2s.z = -Fmenos2s.z;
            Fmenos2s.z *= cosPhi1*icos0;
            //Fmenos2s.z += hij*deta*longitud*icos0;
			//Fmenos2s.z += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*0.5*(u0+u1)*(dy))*longitud*icos0;	//NUEVO
			Fmenos2s.z += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*(u0)*(dy))*longitud*icos0;	//NUEVO
            Fmenos2s.y = u1n*u1*h1*longitud;
            Fmas2s.y = -Fmenos2s.y;
        }
        else {
            Fmenos2s.z = u0n*W0_rot.y*longitud*normal;
            Fmas2s.z = -Fmenos2s.z*cosPhi0*icos1;
            //Fmenos2s.z += hij*deta*longitud*icos0;
			//Fmenos2s.z += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*0.5*(u0+u1)*(dy))*longitud*icos0;	//NUEVO
			Fmenos2s.z += hij*(deta+0.5*(cosPhi0+cosPhi1)*f*(u0)*(dy))*longitud*icos0;	//NUEVO
            Fmenos2s.y = u0n*u0*h0*longitud;
            Fmas2s.y = -Fmenos2s.y;
        }

        cd = 0.5*(1.0 - cvis/(1.0 + exp(-20.0*(1.0-hpos/hij))));
		//cd=0.5;
        fl1.x = W1_rot.y*cosPhi1;
        fl1.y = u1n*fl1.x;
        fl0.x = W0_rot.y*cosPhi0;
        fl0.y = u0n*fl0.x;
        dq = fl1.x-fl0.x;

        des = cd*val*max_autovalor*max_autovalor*dq*longitud;
        Fmenos2s.z -= des*icos0;
        Fmas2s.z += des*icos1;

		if (alfa_cero)
			xiL = xiR = 0.0;
		else {
	        xiL = limitador_monotono(p0m,p0,p1,p1p,aut1);
	        xiR = limitador_monotono(p0m,p0,p1,p1p,aut2);
		}

        if (fabs(aut1-aut2) >= EPSILON) {
            sgL = 1.0*(aut1>0.0) - 1.0*(aut1<0.0);
            sgR = 1.0*(aut2>0.0) - 1.0*(aut2<0.0);
            alr = sgL*(1.0-xiL) + val*aut1*xiL;
            blr = sgR*(1.0-xiR) + val*aut2*xiR;
			aux = 1.0/(aut2-aut1);
            a0 = (alr-blr)*(aut1*aut2)*aux;
            a1 = (blr*aut2-alr*aut1)*aux;
        }
        else {
            a0 = max_autovalor;
            a1 = 0.0;
        }

        Fmenos2->x = Fmas2->x = 0.5*dq;
        //Fmenos2->y = Fmas2->y = 0.5*(fl1.y-fl0.y + hij*deta);
		Fmenos2->y = Fmas2->y = 0.5*(fl1.y-fl0.y + hij*(deta+0.5*(cosPhi0+cosPhi1)*f*0.5*(u0+u1)*(dy)));	//NUEVO
        //DES.x = 0.5*a0*deta;
		DES.x = 0.5*a0*(deta0+0.5*(cosPhi0+cosPhi1)*f*0.5*(u0+u1)*(dy))+a1*Fmenos2->x;	//NUEVO
        DES.y = 0.5*a0*dq + a1*Fmenos2->y;
        Fmenos2->x += fl0.x-DES.x;
        Fmenos2->y += fl0.y-DES.y;
        Fmas2->x += DES.x-fl1.x;
        Fmas2->y += DES.y-fl1.y;
        Fmenos2->x *= longitud;
        Fmenos2->y *= longitud;
        Fmas2->x *= longitud;
        Fmas2->y *= longitud;

        h0m = W0m->x;
        h1p = W1p->x;
        q0mt = -W0m->y*normal;
        q0t = -W0->y*normal;
        q1t = -W1->y*normal;
        q1pt = -W1p->y*normal;
		maxheps = max(h0m,epsilon_h);
        p0m = M_SQRT2*h0m*q0mt / sqrt(h0m*h0m*h0m*h0m + maxheps*maxheps*maxheps*maxheps);
        p0 = u0n = h0rg*q0t;
        p1 = u1n = h1rg*q1t;
		maxheps = max(h1p,epsilon_h);
        p1p = M_SQRT2*h1p*q1pt / sqrt(h1p*h1p*h1p*h1p + maxheps*maxheps*maxheps*maxheps);

        a0 = -1.0*(Fmenos2->x<0.0) + 1.0*(Fmenos2->x>0.0);
        uijn = fabs(uijn)*a0;
		if (alfa_cero)
			xiL = 0.0;
		else
	        xiL = limitador_monotono(p0m, p0, p1, p1p, uijn);
        uijt = 0.5*(u0n*icos0 + u1n*icos1) - 0.5*((1.0-xiL)*a0 + xiL*uijn*val)*(u1n*icos1 - u0n*icos0);
        aux = Fmenos2->x*uijt;
        Fmas3.y = aux*normal;
        Fmas3.z = Fmas2->y*normal*icos1;
        Fmenos3.y = -aux*normal;
        Fmenos3.z = Fmenos2->y*normal*icos0;

        if (max_autovalor < epsilon_h)
            max_autovalor += epsilon_h;

		if (alfa_cero)
			alfa = 0.0;
		else
			alfa = (1.0/(1.0 + exp(-15.0*(1.0 - 0.5*hpos/(hij+EPSILON)))))*(hij < hpos) + 1.0*(hij >= hpos);

		Fmenos2->x = alfa*Fmenos2s.y + (1.0-alfa)*Fmenos3.y;
		Fmenos2->y = alfa*Fmenos2s.z + (1.0-alfa)*Fmenos3.z;

		Fmas2->x = alfa*Fmas2s.y + (1.0-alfa)*Fmas3.y;
		Fmas2->y = alfa*Fmas2s.z + (1.0-alfa)*Fmas3.z;

		a0 = longitud*max_autovalor;
		if (pos_vol0 >= 0) {
			d_acumulador_1[pos_vol0].y += a0;
			acum = d_acumulador_2[pos_vol0];
			acum.x -= Fmenos2->x;
			acum.y -= Fmenos2->y;
			d_acumulador_2[pos_vol0] = acum;
		}

		if (pos_vol1 != -1) {
			d_acumulador_1[pos_vol1].y += a0;
			acum = d_acumulador_2[pos_vol1];
			acum.x -= Fmas2->x;
			acum.y -= Fmas2->y;
			d_acumulador_2[pos_vol1] = acum;
		}
	}
}

__global__ void procesarAristasHorNivel0Paso1GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double *d_anchoVolumenes, int num_volx, int num_voly,
				double borde1, double borde2, double delta_T, double2 *d_acumulador_1,
				double CFL, double epsilon_h, double hpos, double cvis, int tipo)
{
	double2 datos_vol, datos_vol2;
	double normal;
	double Fmas, Fmenos;
	double area0, area1, longitud;
	double cosPhi0, cosPhi1;
	int pos_x_hebra, pos_y_hebra;
	int pos, pos_vol0, pos_vol1;
	int ppv, ppv0, front;
	double H0m, H0, H1, H1p;
	double h0m, h0, h1, h1p;
	double q0y, q1y, aux;
	bool alfa_cero;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_ARI + threadIdx.x;
	pos_y_hebra = 2*(blockIdx.y*NUM_HEBRASY_ARI + threadIdx.y);
	if (tipo == 2) pos_y_hebra++;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra <= num_voly)) {
		pos_vol0 = (pos_y_hebra-1)*num_volx + pos_x_hebra;
		pos_vol1 = pos_vol0 + num_volx;
		longitud = d_anchoVolumenes[pos_y_hebra];
		front = 0;
		normal = 1.0;

		ppv = pos_y_hebra-1;
		ppv0 = ppv*(ppv >= 0);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		h0 = datos_vol.x;
		H0 = datos_vol.y;
		q0y = d_datosVolumenes_2[pos].y;
		datos_vol2 = d_datosVolumenes_3[ppv0];
		area0 = datos_vol2.x;
		cosPhi0 = datos_vol2.y;
		if (ppv < 0) {
			q0y *= borde1;
			pos_vol0 = pos_x_hebra;
			pos_vol1 = -1;
			front = 1;
			normal = -1.0;
		}
		ppv = pos_y_hebra;
		ppv0 = ppv*(ppv < num_voly) + (num_voly-1)*(ppv >= num_voly);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		h1 = datos_vol.x;
		H1 = datos_vol.y;
		q1y = d_datosVolumenes_2[pos].y;
		datos_vol2 = d_datosVolumenes_3[ppv0];
		area1 = datos_vol2.x;
		cosPhi1 = datos_vol2.y;
		if (ppv >= num_voly) {
			q1y *= borde2;
			pos_vol0 = (pos_y_hebra-1)*num_volx + pos_x_hebra;
			pos_vol1 = -1;
		}

		ppv = pos_y_hebra-2;
		ppv0 = ppv*(ppv >= 0);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		h0m = datos_vol.x;
		H0m = datos_vol.y;

		ppv = pos_y_hebra+1;
		ppv0 = ppv*(ppv < num_voly) + (num_voly-1)*(ppv >= num_voly);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		h1p = datos_vol.x;
		H1p = datos_vol.y;

		alfa_cero = ( ((pos_y_hebra < 2) || (pos_y_hebra > num_voly-2)) ? true : false);
		//alfa_cero = ( ((pos_y_hebra < 2) || (pos_y_hebra > num_voly-2)|| (cosPhi0<0.09) ||(cosPhi1<0.09)) ? true : false);
		//alfa_cero = true;
	
		if ((h0 > hpos) && (h1 > hpos) && (h0m > hpos) && (h1p > hpos) && (! alfa_cero)) {
			procesarAristaHor_2s_1step(h0, h1, q0y, longitud, area0, area1, delta_T,
				d_acumulador_1, pos_vol0, pos_vol1, cosPhi0, cosPhi1, &Fmas, &Fmenos);
		}
		else {
			if (front == 1) {
				aux = H0m; H0m = H1p; H1p = aux;
				aux = H0;  H0 = H1;   H1 = aux;
				aux = h0m; h0m = h1p; h1p = aux;
				aux = h0;  h0 = h1;   h1 = aux;
				aux = q0y;  q0y = q1y;   q1y = aux;
				aux = area0; area0 = area1; area1 = aux;
				aux = cosPhi0; cosPhi0 = cosPhi1; cosPhi1 = aux;
			}

			procesarAristaHor_2s_waf_1step(h0m, h0, h1, h1p, q0y, q1y, H0m, H0, H1, H1p, normal,
				longitud, area0, area1, delta_T, d_acumulador_1, pos_vol0, pos_vol1, epsilon_h,
				cosPhi0, cosPhi1, hpos, &Fmas, &Fmenos, alfa_cero);
		}
	}
}

__global__ void procesarAristasHorNivel0Paso2GPU(double2 *d_datosVolumenes_1, double2 *d_datosVolumenes_2,
				double2 *d_datosVolumenes_3, double *d_anchoVolumenes, int num_volx, int num_voly,
				double borde1, double borde2, double delta_T, double2 *d_acumulador_1, double2 *d_acumulador_2,
				double CFL, double epsilon_h, double hpos, double cvis, int tipo, double radio_tierra, double *sinPhi, double T)	//NUEVO
{
	double2 datos_vol, datos_vol2;
	double normal;
	TVec2 Fmas2, Fmenos2;
	double area0, area1, longitud;
	double cosPhi0, cosPhi1;
	int pos_x_hebra, pos_y_hebra;
	int pos, pos_vol0, pos_vol1;
	int ppv, ppv0, front;
	TVec3 W0m, W0, W1, W1p, vaux;
	double H0m, H0, H1, H1p, Haux;
	bool alfa_cero;
	double f;	//NUEVO

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_ARI + threadIdx.x;
	pos_y_hebra = 2*(blockIdx.y*NUM_HEBRASY_ARI + threadIdx.y);
	if (tipo == 2) pos_y_hebra++;

	if ((pos_x_hebra < num_volx) && (pos_y_hebra <= num_voly)) {
		pos_vol0 = (pos_y_hebra-1)*num_volx + pos_x_hebra;
		pos_vol1 = pos_vol0 + num_volx;
		longitud = d_anchoVolumenes[pos_y_hebra];
		front = 0;
		normal = 1.0;
		f = 2*EARTH_ANGULAR_VELOCITY*sinPhi[pos_y_hebra]*T;	//NUEVO (multiplico por T para normalizar el parametro de coriolis)

		ppv = pos_y_hebra-1;
		ppv0 = ppv*(ppv >= 0);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		W0.x = datos_vol.x;
		H0 = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W0.y = datos_vol.x;
		W0.z = datos_vol.y;
		datos_vol2 = d_datosVolumenes_3[ppv0];
		area0 = datos_vol2.x;
		cosPhi0 = datos_vol2.y;
		if (ppv < 0) {
			W0.z *= borde1;
			pos_vol0 = pos_x_hebra;
			pos_vol1 = -1;
			front = 1;
			normal = -1.0;
		}

		ppv = pos_y_hebra;
		ppv0 = ppv*(ppv < num_voly) + (num_voly-1)*(ppv >= num_voly);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		W1.x = datos_vol.x;
		H1 = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W1.y = datos_vol.x;
		W1.z = datos_vol.y;
		datos_vol2 = d_datosVolumenes_3[ppv0];
		area1 = datos_vol2.x;
		cosPhi1 = datos_vol2.y;
		if (ppv >= num_voly) {
			W1.z *= borde2;
			pos_vol0 = (pos_y_hebra-1)*num_volx + pos_x_hebra;
			pos_vol1 = -1;
		}

		ppv = pos_y_hebra-2;
		ppv0 = ppv*(ppv >= 0);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		W0m.x = datos_vol.x;
		H0m = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W0m.y = datos_vol.x;
		W0m.z = datos_vol.y;
            
		ppv = pos_y_hebra+1;
		ppv0 = ppv*(ppv < num_voly) + (num_voly-1)*(ppv >= num_voly);
		pos = ppv0*num_volx + pos_x_hebra;
		datos_vol = d_datosVolumenes_1[pos];
		W1p.x = datos_vol.x;
		H1p = datos_vol.y;
		datos_vol = d_datosVolumenes_2[pos];
		W1p.y = datos_vol.x;
		W1p.z = datos_vol.y;

		alfa_cero = ( ((pos_y_hebra < 2) || (pos_y_hebra > num_voly-2)) ? true : false);
		//alfa_cero = ( ((pos_y_hebra < 2) || (pos_y_hebra > num_voly-2)|| (cosPhi0<0.09) ||(cosPhi1<0.09)) ? true : false);
		//alfa_cero = true;
	
		if ((W0.x > hpos) && (W1.x > hpos) && (W0m.x > hpos) && (W1p.x > hpos) && (! alfa_cero)) {
			procesarAristaHor_2s_2step(&W0, &W1, H0, H1, longitud, area0, area1, delta_T, d_acumulador_1,
				d_acumulador_2, pos_vol0, pos_vol1, CFL, epsilon_h, cosPhi0, cosPhi1, cvis, hpos, &Fmas2, &Fmenos2, radio_tierra, f);	//NUEVO (MODIFICADO)
		}
		else {
			if (front == 1) {
				v_copy3(&W0m, &vaux); v_copy3(&W1p, &W0m); v_copy3(&vaux, &W1p);
				v_copy3(&W0, &vaux);  v_copy3(&W1, &W0);   v_copy3(&vaux, &W1);
				Haux = H0m; H0m = H1p; H1p = Haux;
				Haux = H0;  H0 = H1;   H1 = Haux;
				Haux = area0; area0 = area1; area1 = Haux;
				Haux = cosPhi0; cosPhi0 = cosPhi1; cosPhi1 = Haux;
			}

			procesarAristaHor_2s_waf_2step(&W0m, &W0, &W1, &W1p, H0m, H0, H1, H1p, normal, longitud,
				area0, area1, delta_T, d_acumulador_1, d_acumulador_2, pos_vol0, pos_vol1, CFL,
				epsilon_h, cosPhi0, cosPhi1, cvis, hpos, &Fmas2, &Fmenos2, alfa_cero, radio_tierra, f);	//NUEVO (MODIFICADO)
		}
	}
}

// DELTA_T INICIAL

__device__ void procesarAristaDeltaTInicial(TVec3 *W0, TVec3 *W1, double H0, double H1, double normal1_x, double normal1_y,
			double longitud, double2 *d_acumulador_2, int pos_vol0, int pos_vol1, double epsilon_h, double cosPhi0)
{
	double h1ij, u1ij_n;
	double max_autovalor;
	double u0n, u1n, maxheps;
	double h0, h1, sqrt_h0, sqrt_h1;
	TVec2 W0_rot, W1_rot;
	double2 acum0, acum1;

	W0_rot.x = W0->x;
	W0_rot.y = W0->y*normal1_x + W0->z*normal1_y;
	W1_rot.x = W1->x;
	W1_rot.y = W1->y*normal1_x + W1->z*normal1_y;

	h0 = W0_rot.x;
	h1 = W1_rot.x;
	h1ij = 0.5*(h0 + h1);

	maxheps = max(h0,epsilon_h);
	u0n = M_SQRT2*h0*W0_rot.y / sqrt(h0*h0*h0*h0 + maxheps*maxheps*maxheps*maxheps);
	maxheps = max(h1,epsilon_h);
	u1n = M_SQRT2*h1*W1_rot.y / sqrt(h1*h1*h1*h1 + maxheps*maxheps*maxheps*maxheps);
	sqrt_h0 = sqrt(h0);
	sqrt_h1 = sqrt(h1);
	u1ij_n = (sqrt_h0*u0n + sqrt_h1*u1n) / (sqrt_h0 + sqrt_h1 + EPSILON);

	max_autovalor = fabs(u1ij_n)+sqrt(h1ij);

	if (max_autovalor < epsilon_h)
		max_autovalor += epsilon_h;

	h0 = longitud*max_autovalor/cosPhi0;
	acum0 = d_acumulador_2[pos_vol0];
	acum0.y += h0;
	d_acumulador_2[pos_vol0] = acum0;

	if (pos_vol1 != -1) {
		acum1 = d_acumulador_2[pos_vol1];
		acum1.y += h0;
		d_acumulador_2[pos_vol1] = acum1;
	}
}

__global__ void procesarAristasVerDeltaTInicialNivel0GPU(double2 *d_datosVolumenesNivel0_1, double2 *d_datosVolumenesNivel0_2,
				double2 *d_datosVolumenes_3, int numVolxNivel0, int numVolyNivel0, double borde1, double borde2,
				double *d_altoVolumenesNivel0, double2 *d_acumulador_1, double epsilon_h, int tipo)
{
	double2 datos_vol0, datos_vol1;
	double longitud, cosPhi0;
	double H0, H1;
	int pos_x_hebra, pos_y_hebra;
	int pos, pos_vol0;
	TVec3 W0, W1;

	pos_x_hebra = 2*(blockIdx.x*NUM_HEBRASX_ARI + threadIdx.x);
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_ARI + threadIdx.y;
	if (tipo == 2) pos_x_hebra++;

	if ((pos_x_hebra <= numVolxNivel0) && (pos_y_hebra < numVolyNivel0)) {
		longitud = d_altoVolumenesNivel0[pos_y_hebra];
		cosPhi0 = d_datosVolumenes_3[pos_y_hebra].y;
		if (pos_x_hebra == 0) {
			pos = pos_y_hebra*numVolxNivel0 + pos_x_hebra;
			datos_vol0 = d_datosVolumenesNivel0_1[pos];
			W0.x = datos_vol0.x;
			H0 = datos_vol0.y;
			datos_vol0 = d_datosVolumenesNivel0_2[pos];
			W0.y = datos_vol0.x;
			W0.z = datos_vol0.y;
			//NUEVO
			if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
				pos = pos_y_hebra*numVolxNivel0 + numVolxNivel0 - 1;
				datos_vol1 = d_datosVolumenesNivel0_1[pos];
				W1.x = datos_vol1.x;
				H1 = datos_vol1.y;
				datos_vol1 = d_datosVolumenesNivel0_2[pos];
				W1.y = datos_vol1.x;
				W1.z = datos_vol1.y;
			}
			else{
				W1.x = W0.x;
				W1.y = W0.y*borde1;
				W1.z = W0.z;
				H1 = H0;
			}
			pos_vol0 = pos_y_hebra*numVolxNivel0;
			procesarAristaDeltaTInicial(&W0, &W1, H0, H1, -1.0, 0.0, longitud,
				d_acumulador_1, pos_vol0, -1, epsilon_h, cosPhi0);
		}
		else {
			pos = pos_y_hebra*numVolxNivel0 + pos_x_hebra-1;
			datos_vol0 = d_datosVolumenesNivel0_1[pos];
			W0.x = datos_vol0.x;
			H0 = datos_vol0.y;
			datos_vol0 = d_datosVolumenesNivel0_2[pos];
			W0.y = datos_vol0.x;
			W0.z = datos_vol0.y;
			if (pos_x_hebra == numVolxNivel0) {
				//NUEVO
				if (fabs(borde1-2.0)<EPSILON && fabs(borde2-2.0)<EPSILON){
					pos = pos_y_hebra*numVolxNivel0;
					datos_vol1 = d_datosVolumenesNivel0_1[pos];
					W1.x = datos_vol1.x;
					H1 = datos_vol1.y;
					datos_vol1 = d_datosVolumenesNivel0_2[pos];
					W1.y = datos_vol1.x;
					W1.z = datos_vol1.y;
				}
				else{
					W1.x = W0.x;
					W1.y = W0.y*borde2;
					W1.z = W0.z;
					H1 = H0;
				}
				pos_vol0 = (pos_y_hebra+1)*numVolxNivel0 - 1;
				procesarAristaDeltaTInicial(&W0, &W1, H0, H1, 1.0, 0.0, longitud,
					d_acumulador_1, pos_vol0, -1, epsilon_h, cosPhi0);
			}
			else {
				pos = pos_y_hebra*numVolxNivel0 + pos_x_hebra;
				datos_vol1 = d_datosVolumenesNivel0_1[pos];
				W1.x = datos_vol1.x;
				H1 = datos_vol1.y;
				datos_vol1 = d_datosVolumenesNivel0_2[pos];
				W1.y = datos_vol1.x;
				W1.z = datos_vol1.y;

				pos_vol0 = pos_y_hebra*numVolxNivel0 + pos_x_hebra - 1;
				procesarAristaDeltaTInicial(&W0, &W1, H0, H1, 1.0, 0.0, longitud,
					d_acumulador_1, pos_vol0, pos_vol0+1, epsilon_h, cosPhi0);
			}
		}
	}
}

__global__ void procesarAristasHorDeltaTInicialNivel0GPU(double2 *d_datosVolumenesNivel0_1, double2 *d_datosVolumenesNivel0_2,
				int numVolxNivel0, int numVolyNivel0, double borde1, double borde2, double *d_anchoVolumenesNivel0,
				double2 *d_acumulador_1, double epsilon_h, int tipo)
{
	double2 datos_vol0, datos_vol1;
	double H0, H1;
	double longitud;
	int pos_x_hebra, pos_y_hebra;
	int pos, pos_vol0, pos_vol1;
	TVec3 W0, W1;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_ARI + threadIdx.x;
	pos_y_hebra = 2*(blockIdx.y*NUM_HEBRASY_ARI + threadIdx.y);
	if (tipo == 2) pos_y_hebra++;

	if ((pos_x_hebra < numVolxNivel0) && (pos_y_hebra <= numVolyNivel0)) {
		longitud = d_anchoVolumenesNivel0[pos_y_hebra];
		if (pos_y_hebra == 0) {
			pos = pos_y_hebra*numVolxNivel0 + pos_x_hebra;
			datos_vol0 = d_datosVolumenesNivel0_1[pos];
			W0.x = datos_vol0.x;
			H0 = datos_vol0.y;
			datos_vol0 = d_datosVolumenesNivel0_2[pos];
			W0.y = datos_vol0.x;
			W0.z = datos_vol0.y;
			W1.x = W0.x;
			W1.y = W0.y;
			W1.z = W0.z*borde1;

			procesarAristaDeltaTInicial(&W0, &W1, H0, H0, 0.0, -1.0, longitud,
				d_acumulador_1, pos_x_hebra, -1, epsilon_h, 1.0);
		}
		else {
			pos = (pos_y_hebra-1)*numVolxNivel0 + pos_x_hebra;
			datos_vol0 = d_datosVolumenesNivel0_1[pos];
			W0.x = datos_vol0.x;
			H0 = datos_vol0.y;
			datos_vol0 = d_datosVolumenesNivel0_2[pos];
			W0.y = datos_vol0.x;
			W0.z = datos_vol0.y;
			if (pos_y_hebra == numVolyNivel0) {
				W1.x = W0.x;
				W1.y = W0.y;
				W1.z = W0.z*borde2;

				pos_vol0 = (pos_y_hebra-1)*numVolxNivel0 + pos_x_hebra;
				procesarAristaDeltaTInicial(&W0, &W1, H0, H0, 0.0, 1.0, longitud,
					d_acumulador_1, pos_vol0, -1, epsilon_h, 1.0);
			}
			else {
				pos = pos_y_hebra*numVolxNivel0 + pos_x_hebra;
				datos_vol1 = d_datosVolumenesNivel0_1[pos];
				W1.x = datos_vol1.x;
				H1 = datos_vol1.y;
				datos_vol1 = d_datosVolumenesNivel0_2[pos];
				W1.y = datos_vol1.x;
				W1.z = datos_vol1.y;

				pos_vol1 = pos_y_hebra*numVolxNivel0 + pos_x_hebra;
				pos_vol0 = pos_vol1 - numVolxNivel0;
				procesarAristaDeltaTInicial(&W0, &W1, H0, H1, 0.0, 1.0, longitud,
					d_acumulador_1, pos_vol0, pos_vol1, epsilon_h, 1.0);
			}
		}
	}
}

#endif
