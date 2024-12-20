#ifndef _DEFORMACION_DINAMICA_H_
#define _DEFORMACION_DINAMICA_H_

#include <stdio.h>
#include "Constantes.hxx"
#include "Deformacion.cu"
#include "netcdf.cu"

// NetCDF external functions
//extern int abrirGRDDefDinamica(const char *nombre_fich);
//extern void leerEstadoDefDinamicaGRD(int nvx, int nvy, int num, float *def);
//extern void cerrarGRDDefDinamica();

//__global__ void interpolarDefEnTiempoNivel0GPU(double *d_defSigNivel0, double defTime_ant, double defTime_sig,
//			double tiempoAnt, double tiempoAct, double *d_deformacionInterpoladaNivel0, int nvxDef, int nvyDef)
//{
//	int pos_x_hebra, pos_y_hebra;
//	int pos_def;
//	double peso, U_Z;
//
//	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
//	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;
//
//	if ((pos_x_hebra < nvxDef) && (pos_y_hebra < nvyDef)) {
//		pos_def = pos_y_hebra*nvxDef + pos_x_hebra;
//		U_Z = d_defSigNivel0[pos_def];
//
//		if (fabs(defTime_sig-defTime_ant) > EPSILON) {
//			peso = (tiempoAct - tiempoAnt) / (defTime_sig - defTime_ant);
//			U_Z = peso*U_Z;
//		}
//
//		d_deformacionInterpoladaNivel0[pos_def] = U_Z;
//	}
//}

void normalizaPresion(float *presion0, int nvxDef, int nvyDef, double H, double Q)	//NUEVO
{
	int i, j, pos;
	double normalized_pressure;
	double U = Q/H;	//NUEVO

	for (j=0; j<nvyDef; j++) {
		for (i=0; i<nvxDef; i++) {
			pos = j*nvxDef + i;
			normalized_pressure = (presion0[pos])/(U*U*WATER_DENSITY); //NUEVO
			presion0[pos] = (float) normalized_pressure; //aqui guardo el campo de presiones inicial
		}
	}
}


//NUEVO
__global__ void copiarBathyDesdeComponenteY(double *d_original_bathy, double2 *d_datosVolumenesNivel_1, int nvx, int nvy) {
	int pos_x_hebra, pos_y_hebra;
	int pos;

    pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

    if ((pos_x_hebra < nvx) && (pos_y_hebra < nvy)) {
		pos = pos_y_hebra*nvx + pos_x_hebra;
        d_original_bathy[pos] = d_datosVolumenesNivel_1[pos].y;
    }
}

//NUEVO
__global__ void copiarBathyHaciaComponenteY(double *d_original_bathy, double2 *d_datosVolumenesNivel_1, int nvx, int nvy) {
	int pos_x_hebra, pos_y_hebra;
	int pos;

    pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

    if ((pos_x_hebra < nvx) && (pos_y_hebra < nvy)) {
		pos = pos_y_hebra*nvx + pos_x_hebra;
        d_datosVolumenesNivel_1[pos].y = d_original_bathy[pos];
    }
}

__global__ void copiarPresion1HaciaPresion0(float *d_presion0, float *d_presion1, int nvx, int nvy) {
	int pos_x_hebra, pos_y_hebra;
	int pos;

    pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

    if ((pos_x_hebra < nvx) && (pos_y_hebra < nvy)) {
		pos = pos_y_hebra*nvx + pos_x_hebra;
        d_presion0[pos] = d_presion1[pos];
    }
}


__global__ void interpolarDefEnTiempoNivel0GPU(float *d_defAntNivel0, float *d_defSigNivel0, double defTime_ant, double defTime_sig,
			double tiempoAnt, double tiempoAct, double *d_deformacionInterpoladaNivel0, int nvxDef, int nvyDef) //NUEVO
{
	int pos_x_hebra, pos_y_hebra;
	int pos_def;
	double peso, U_Z_sig, U_Z_ant; //NUEVO

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < nvxDef) && (pos_y_hebra < nvyDef)) {
		pos_def = pos_y_hebra*nvxDef + pos_x_hebra;
		U_Z_sig = d_defSigNivel0[pos_def];
		U_Z_ant = d_defAntNivel0[pos_def];

		if (fabs(defTime_sig-defTime_ant) > EPSILON) {
			peso = (defTime_sig-tiempoAct) / (defTime_sig - defTime_ant); //NUEVO
		}
		else{
			peso = 1.0; //aqui a lo mejor hay que poner 0.0
		}
		d_deformacionInterpoladaNivel0[pos_def] = peso*U_Z_ant + (1-peso)*U_Z_sig;
	}
}

__global__ void sumarDeformacionADatosGPU(double2 *d_datosVolumenesNivel0_1, double *d_deformacionNivel0,
				int nvxNivel0, int nvyNivel0)
{
	int pos_x_hebra, pos_y_hebra;
	int pos_datos;
	double U_Z;

	pos_x_hebra = blockIdx.x*NUM_HEBRASX_EST + threadIdx.x;
	pos_y_hebra = blockIdx.y*NUM_HEBRASY_EST + threadIdx.y;

	if ((pos_x_hebra < nvxNivel0) && (pos_y_hebra < nvyNivel0)) {
		pos_datos = pos_y_hebra*nvxNivel0 + pos_x_hebra;
		U_Z = d_deformacionNivel0[pos_datos];
		d_datosVolumenesNivel0_1[pos_datos].y -= U_Z;
	}
}

//ASI ESTABA ANTES
//void aplicarDefDinamica(int indiceFallaSig, int nvxNivel0, int nvyNivel0, double2 *d_datosVolumenesNivel_1,
//		double **d_deformacionNivel0, double *d_deformacionInterpoladaNivel0, double *defTime,
//		double tiempoAnt, double tiempoAct, dim3 blockGridEstNivel, dim3 threadBlockEst)
//{
//	double *d_defSigNivel0, *d_defAntNivel0; //NUEVO
//	double defTime_ant, defTime_sig;
//	int indiceFallaAnt = max(0, indiceFallaSig-1);
//
//	d_defAntNivel0 = d_deformacionNivel0[indiceFallaAnt]; //NUEVO
//	d_defSigNivel0 = d_deformacionNivel0[indiceFallaSig];
//	defTime_ant = defTime[indiceFallaAnt];
//	defTime_sig = defTime[indiceFallaSig];
//	interpolarDefEnTiempoNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_defAntNivel0, d_defSigNivel0, defTime_ant,
//		defTime_sig, tiempoAnt, tiempoAct, d_deformacionInterpoladaNivel0, nvxNivel0, nvyNivel0); //NUEVO
//
//	sumarDeformacionADatosGPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1,
//		d_deformacionInterpoladaNivel0, nvxNivel0, nvyNivel0);
//}

void aplicarDefDinamica(int indiceFallaSig, int nvxNivel0, int nvyNivel0, double2 *d_datosVolumenesNivel_1,
		float *d_presion0, float *d_presion1, double *d_deformacionInterpoladaNivel0, double *defTime,
		double tiempoAnt, double tiempoAct, dim3 blockGridEstNivel, dim3 threadBlockEst)
{
	float *d_defSigNivel0, *d_defAntNivel0; //NUEVO
	double defTime_ant, defTime_sig;
	int indiceFallaAnt = max(0, indiceFallaSig-1);


	d_defAntNivel0 = d_presion0;
	d_defSigNivel0 = d_presion1;

	//d_defAntNivel0 = d_deformacionNivel0[indiceFallaAnt]; //NUEVO
	//d_defSigNivel0 = d_deformacionNivel0[indiceFallaSig];
	defTime_ant = defTime[indiceFallaAnt];
	defTime_sig = defTime[indiceFallaSig];
	interpolarDefEnTiempoNivel0GPU<<<blockGridEstNivel, threadBlockEst>>>(d_defAntNivel0, d_defSigNivel0, defTime_ant,
		defTime_sig, tiempoAnt, tiempoAct, d_deformacionInterpoladaNivel0, nvxNivel0, nvyNivel0); //NUEVO

	sumarDeformacionADatosGPU<<<blockGridEstNivel, threadBlockEst>>>(d_datosVolumenesNivel_1,
		d_deformacionInterpoladaNivel0, nvxNivel0, nvyNivel0);
}

void comprobarYAplicarDeformacionDinamica(double2 *d_datosVolumenesNivel_1, float *d_presion0, float *d_presion1,
		double *d_deltaTVolumenesNivel, double tiempoAntSubmalla, double tiempoActSubmalla, int nvxNivel0,
		int nvyNivel0, int okada_flag, int *fallaOkada, double *defTime, int numEstadosDefDinamica,
		int *indiceEstadoSigDefDim, dim3 blockGridEstNivel, dim3 threadBlockEst, double T, float *presion0, string directorio, string fich_def, double H, double Q)
{
	if ((*fallaOkada == 0) && (okada_flag == DYNAMIC_DEFORMATION)) {
		if ((tiempoActSubmalla*T >= defTime[0]) && (tiempoActSubmalla*T <= defTime[numEstadosDefDinamica-1])) {
			fprintf(stdout, "Applying pressure forcing\n");
			aplicarDefDinamica(*indiceEstadoSigDefDim, nvxNivel0, nvyNivel0, d_datosVolumenesNivel_1, d_presion0, d_presion1,
				d_deltaTVolumenesNivel, defTime, tiempoAntSubmalla*T, tiempoActSubmalla*T, blockGridEstNivel, threadBlockEst);
			if (tiempoActSubmalla*T >= defTime[*indiceEstadoSigDefDim]) {
				*indiceEstadoSigDefDim = min((*indiceEstadoSigDefDim)+1, numEstadosDefDinamica-1);
				copiarPresion1HaciaPresion0<<<blockGridEstNivel, threadBlockEst>>>(d_presion0, d_presion1, nvxNivel0, nvyNivel0);
				abrirGRDDefDinamica((directorio+fich_def).c_str());
				leerEstadoDefDinamicaGRD(nvxNivel0, nvyNivel0, *indiceEstadoSigDefDim, presion0);	//en presion0 se guarda el siguiente instante de presion en CPU, que se copiara a GPU
				cerrarGRDDefDinamica();
				normalizaPresion(presion0, nvxNivel0, nvyNivel0, H, Q);
				int64_t tam2 = (int64_t) nvxNivel0*nvyNivel0*sizeof(float);
				cudaMemcpy(d_presion1, presion0, tam2, cudaMemcpyHostToDevice);
			}
			if (tiempoActSubmalla*T >= defTime[numEstadosDefDinamica-1]) {
				(*fallaOkada)++;
			}
		}
	}
}




#endif
