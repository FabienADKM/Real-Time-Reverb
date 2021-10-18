/******************************************/
/*
  duplex.cpp
  by Gary P. Scavone, 2006-2007.

  This program opens a duplex stream and passes
  input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
//#include <sys/time.h>
//#include <sys/resource.h>
#include <math.h>
////////
#include <windows.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define PI 3.14159

double* _sintbl = 0;
int maxfftsize = 0;
int fft(double* x, double* y, const int m);
int ifft(double* x, double* y, const int m);
int fftr(double* x, double* y, const int m);
int ifftr(double* x, double* y, const int l);
static int checkm(const int m);
int get_nextpow2(int n);
char* getmem(int leng, unsigned size);
double* dgetmem(int leng);







///////////////////////////////
// FFT functions
int get_nextpow2(int n)
{
    int k = 1;
    while (k < n) {
        k *= 2;
    }

    return k;
}

int fftr(double* x, double* y, const int m)
{
    int i, j;
    double* xp, * yp, * xq;
    double* yq;
    int mv2, n, tblsize;
    double xt, yt, * sinp, * cosp;
    double arg;

    mv2 = m / 2;

    /* separate even and odd  */
    xq = xp = x;
    yp = y;
    for (i = mv2; --i >= 0;) {
        *xp++ = *xq++;
        *yp++ = *xq++;
    }

    if (fft(x, y, mv2) == -1)    /* m / 2 point fft */
        return (-1);


    /***********************
    * SIN table generation *
    ***********************/

    if ((_sintbl == 0) || (maxfftsize < m)) {
        tblsize = m - m / 4 + 1;
        arg = PI / m * 2;
        if (_sintbl != 0)
            free(_sintbl);
        _sintbl = sinp = dgetmem(tblsize);
        *sinp++ = 0;
        for (j = 1; j < tblsize; j++)
            *sinp++ = sin(arg * (double)j);
        _sintbl[m / 2] = 0;
        maxfftsize = m;
    }
    //printf("Debug: m=%i, maxfftsize=%i\n",m,maxfftsize);

    n = maxfftsize / m;
    sinp = _sintbl;
    cosp = _sintbl + maxfftsize / 4;

    xp = x;
    yp = y;
    xq = xp + m;
    yq = yp + m;
    *(xp + mv2) = *xp - *yp;
    *xp = *xp + *yp;
    *(yp + mv2) = *yp = 0;

    for (i = mv2, j = mv2 - 2; --i; j -= 2) {
        ++xp;
        ++yp;
        sinp += n;
        cosp += n;
        yt = *yp + *(yp + j);
        xt = *xp - *(xp + j);
        *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
        *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
    }

    xp = x + 1;
    yp = y + 1;
    xq = x + m;
    yq = y + m;

    for (i = mv2; --i;) {
        *xp++ = *(--xq);
        *yp++ = -(*(--yq));
    }

    return (0);
}

int ifftr(double* x, double* y, const int l)
{
    int i;
    double* xp, * yp;

    fftr(x, y, l);

    xp = x;
    yp = y;
    i = l;
    while (i--) {
        *xp++ /= l;
        *yp++ /= -l;
    }

    return (0);
}


static int checkm(const int m)
{
    int k;

    for (k = 4; k <= m; k <<= 1) {
        if (k == m)
            return (0);
    }
    fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n", m);

    return (-1);
}

int fft(double* x, double* y, const int m)
{
    int j, lmx, li;
    double* xp, * yp;
    double* sinp, * cosp;
    int lf, lix, tblsize;
    int mv2, mm1;
    double t1, t2;
    double arg;
    int checkm(const int);

    /**************
    * RADIX-2 FFT *
    **************/

    if (checkm(m))
        return (-1);

    /***********************
    * SIN table generation *
    ***********************/

    if ((_sintbl == 0) || (maxfftsize < m)) {
        tblsize = m - m / 4 + 1;
        arg = PI / m * 2;
        if (_sintbl != 0)
            free(_sintbl);
        _sintbl = sinp = dgetmem(tblsize);
        *sinp++ = 0;
        for (j = 1; j < tblsize; j++)
            *sinp++ = sin(arg * (double)j);
        _sintbl[m / 2] = 0;
        maxfftsize = m;
    }

    lf = maxfftsize / m;
    lmx = m;

    for (;;) {
        lix = lmx;
        lmx /= 2;
        if (lmx <= 1)
            break;
        sinp = _sintbl;
        cosp = _sintbl + maxfftsize / 4;
        for (j = 0; j < lmx; j++) {
            xp = &x[j];
            yp = &y[j];
            for (li = lix; li <= m; li += lix) {
                t1 = *(xp)-*(xp + lmx);
                t2 = *(yp)-*(yp + lmx);
                *(xp) += *(xp + lmx);
                *(yp) += *(yp + lmx);
                *(xp + lmx) = *cosp * t1 + *sinp * t2;
                *(yp + lmx) = *cosp * t2 - *sinp * t1;
                xp += lix;
                yp += lix;
            }
            sinp += lf;
            cosp += lf;
        }
        lf += lf;
    }

    xp = x;
    yp = y;
    for (li = m / 2; li--; xp += 2, yp += 2) {
        t1 = *(xp)-*(xp + 1);
        t2 = *(yp)-*(yp + 1);
        *(xp) += *(xp + 1);
        *(yp) += *(yp + 1);
        *(xp + 1) = t1;
        *(yp + 1) = t2;
    }

    /***************
    * bit reversal *
    ***************/
    j = 0;
    xp = x;
    yp = y;
    mv2 = m / 2;
    mm1 = m - 1;
    for (lmx = 0; lmx < mm1; lmx++) {
        if ((li = lmx - j) < 0) {
            t1 = *(xp);
            t2 = *(yp);
            *(xp) = *(xp + li);
            *(yp) = *(yp + li);
            *(xp + li) = t1;
            *(yp + li) = t2;
        }
        li = mv2;
        while (li <= j) {
            j -= li;
            li /= 2;
        }
        j += li;
        xp = x + j;
        yp = y + j;
    }

    return (0);
}

int ifft(double* x, double* y, const int m)
{
    int i;

    if (fft(y, x, m) == -1)
        return (-1);

    for (i = m; --i >= 0; ++x, ++y) {
        *x /= m;
        *y /= m;
    }

    return (0);
}

double* dgetmem(int leng)
{
    return ((double*)getmem(leng, sizeof(double)));
}

char* getmem(int leng, unsigned size)
{
    char* p = NULL;

    if ((p = (char*)calloc(leng, size)) == NULL) {
        fprintf(stderr, "Memory allocation error !\n");
        exit(3);
    }
    return (p);
}


/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8
*/

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

typedef struct {
    int bufferBytes;
    unsigned int bufferFrames;
    int hSize;
    int Ly;
    int Ly2;
    MY_TYPE *h;
    MY_TYPE *temp;
    MY_TYPE *conv;
    MY_TYPE *x_fft_result_re;
    MY_TYPE *x_fft_result_im;
    MY_TYPE *h_fft_result_re;
    MY_TYPE *h_fft_result_im;
    MY_TYPE *outputbuffer;
    
}paData;

/*
typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64
*/

void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device to use (default = 0),\n";
  std::cout << "    oDevice = optional output device to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
  exit( 0 );
}


int time_conv(MY_TYPE* inputBuffer, MY_TYPE* h, MY_TYPE* conv,  int hSize,  int BufferFrames) {
    int kmin = 0;
    int kmax = 0;
    unsigned int L = BufferFrames;  // taille du signal x
    unsigned int M = hSize;
    int tmp = 0;
    MY_TYPE* x = inputBuffer;
    for (int n = 0; n <L + M-1; n++) {
        // on fixe kmin et kmax
        tmp = 0;
        kmin = MAX(int(n - M + 1), 0);
        kmax = MIN(n, L);

        // on calcule la convolution pour chaque échantillon
        for (int k = kmin; k < kmax; k++) {
            // tmp = tmp + *((int *)data+k)* (*((int *)data+n-k+1));
            tmp = tmp + x[k] * h[n - k + 1];
        }
        conv[n] = tmp;
    }

    return 0;
}

int freq_conv(MY_TYPE* conv, MY_TYPE* x_fftresult_r, MY_TYPE* x_fftresult_im, MY_TYPE* h_fftresult_r,
    MY_TYPE* h_fftresult_im, int Ly2, int Ly) {

    fftr(x_fftresult_r, x_fftresult_im, Ly2);
    for (int n = 0; n < Ly2; n++) {
    x_fftresult_r[n] = x_fftresult_r[n] * h_fftresult_r[n] - x_fftresult_im[n] * h_fftresult_im[n];
    x_fftresult_im[n] = x_fftresult_r[n] * h_fftresult_im[n] + h_fftresult_r[n]* x_fftresult_im[n];
    }
    ifft(x_fftresult_r, x_fftresult_im, Ly2);
    memcpy(conv, x_fftresult_r, sizeof(MY_TYPE) * Ly);
    //Réinitialisation des tableaux à zéro
    for (int n = 0; n < Ly2; n++) {
        x_fftresult_r[n] = 0;
        x_fftresult_im[n] = 0;
    }
    return 0;
}



int inout( void *outputBuffer, void *inputBuffer, unsigned int /*nBufferFrames*/,
           double /*streamTime*/, RtAudioStreamStatus status, void* userData )
{  
  //Mesure du temps d'execution 
  long depart = 0;
  long fin = 0;
  long duree = 0;
  depart = GetTickCount();

  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;
  paData* p_userdata = (paData*)userData; 
  //Traitement pour la convolution fréquentielle
  memcpy(p_userdata->x_fft_result_re, inputBuffer, p_userdata->bufferBytes);
  
  //Convolution temporelle
  //if(time_conv((MY_TYPE *)inputBuffer, p_userdata->h, p_userdata->conv, p_userdata->hSize, p_userdata->bufferFrames) != 0){ fputs("Convolution error", stderr); exit(7); }
  //Convolution fréquentielle
  if(freq_conv(p_userdata->conv, p_userdata->x_fft_result_re, p_userdata->x_fft_result_im, p_userdata->h_fft_result_re, p_userdata->h_fft_result_im, p_userdata->Ly2, p_userdata->Ly) !=0){ fputs("Convolution error", stderr); exit(7); }
  
  for (int i = 0; i < p_userdata->Ly - p_userdata->bufferFrames - 1; i++) {
      p_userdata->conv[i] = p_userdata->conv[i] + p_userdata->temp[i+p_userdata->bufferFrames];
  }
  for (int i = 0; i < p_userdata->bufferFrames; ++i) {
      p_userdata->outputbuffer[i] = p_userdata->conv[i];
  }
  for (int i = 0; i < p_userdata->Ly - 1; i++) {
      p_userdata->temp[i] = p_userdata->conv[i];
  }

  memcpy( outputBuffer, p_userdata->outputbuffer, p_userdata->bufferBytes);
  
  // Code à mesurer
  fin = GetTickCount();
  duree = fin - depart;
  printf("time: %ld ms\n", duree);
  return 0;
}

int main( int argc, char *argv[] )
{
  double* _sintbl = 0;
  int maxfftsize = 0;
  unsigned int channels, fs, bufferBytes, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

  // Minimal command-line checking
  if (argc < 3 || argc > 7 ) usage();

  RtAudio adac;
  if ( adac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  if ( argc > 3 )
    iDevice = (unsigned int) atoi(argv[3]);
  if ( argc > 4 )
    oDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    iOffset = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    oOffset = (unsigned int) atoi(argv[6]);

  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 800;
  RtAudio::StreamParameters iParams, oParams;
  iParams.deviceId = iDevice;
  iParams.nChannels = channels;
  iParams.firstChannel = iOffset;
  oParams.deviceId = oDevice;
  oParams.nChannels = channels;
  oParams.firstChannel = oOffset;
  
  if ( iDevice == 0 )
    iParams.deviceId = adac.getDefaultInputDevice();
  if ( oDevice == 0 )
    oParams.deviceId = adac.getDefaultOutputDevice();

  RtAudio::StreamOptions options;
  //options.flags |= RTAUDIO_NONINTERLEAVED;

  bufferBytes = bufferFrames * channels * sizeof( MY_TYPE );

  /*Open impulsionnel response*/
  FILE* pFile;
  int hSize;
  size_t result;

  pFile = fopen("impres", "rb");
  if (pFile == NULL) { fputs("File error", stderr); exit(1); }

  /*
  // obtain impulse size:
  fseek(pFile, 0, SEEK_END);
  hSize = ftell(pFile)/8;
  std::cout << "\n Valeur de hSize = " << hSize << ").\n";
  rewind(pFile);
  */
  hSize = 20000;

  /*Declaration of userData*/
  paData* userData = NULL;
  userData = (paData*)malloc(sizeof(paData));

  /*Initialisation of userData*/
  int Ly = bufferFrames + hSize - 1;
  int Ly2 = get_nextpow2(Ly);

  if (userData != NULL) {
      userData->bufferFrames = bufferFrames;
      userData->bufferBytes = bufferBytes;
      userData->hSize = hSize;
      userData->Ly = Ly;
      userData->Ly2 = Ly2;
      userData->h = (MY_TYPE*)malloc(hSize * sizeof(MY_TYPE));
      if (userData->h == NULL) { fputs("Memory error", stderr); exit(2); }
      result = fread(userData->h, sizeof(MY_TYPE), hSize, pFile);
      if (result != hSize) { fputs("Reading error", stderr); exit(3); }
      userData->temp = (MY_TYPE*)calloc(Ly, sizeof(MY_TYPE));
      if (userData->temp == NULL) { fputs("Memory error", stderr); exit(4); }
      userData->conv = (MY_TYPE*)calloc(Ly, sizeof(MY_TYPE));
      if (userData->conv == NULL) { fputs("Memory error", stderr); exit(5); }
      userData->x_fft_result_re = (MY_TYPE*)calloc(Ly2, sizeof(MY_TYPE));
      if (userData->x_fft_result_re == NULL) { fputs("Memory error", stderr); exit(5); }
      userData->x_fft_result_im = (MY_TYPE*)calloc(Ly2, sizeof(MY_TYPE));
      if (userData->x_fft_result_im == NULL) { fputs("Memory error", stderr); exit(5); }
      userData->h_fft_result_re = (MY_TYPE*)calloc(Ly2, sizeof(MY_TYPE));
      if (userData->h_fft_result_re == NULL) { fputs("Memory error", stderr); exit(5); }
      memcpy(userData->h_fft_result_re, userData->h, hSize * sizeof(MY_TYPE));
      userData->h_fft_result_im = (MY_TYPE*)calloc(Ly2, sizeof(MY_TYPE));
      if (userData->h_fft_result_im == NULL) { fputs("Memory error", stderr); exit(5); }
      
      

      //calcul de la fft de h
      if(fftr(userData->h_fft_result_re , userData->h_fft_result_im , Ly2) != 0) { fputs("fftr fail", stderr); exit(5); };
      
      userData->outputbuffer = (MY_TYPE*)calloc(bufferFrames, sizeof(MY_TYPE));
      if (userData->outputbuffer == NULL) { fputs("Memory error", stderr); exit(5); }
  }
  // closing file
  fclose(pFile);
  
  


  try {
    adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, (void*)userData, &options );
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    exit( 1 );
  }

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  try {
    adac.startStream();

    char input;
    std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
    std::cin.get(input);

    // Stop the stream.
    adac.stopStream();
  }
  catch ( RtAudioError& e ) {
    std::cout << '\n' << e.getMessage() << '\n' << std::endl;
    goto cleanup;
  }
  free(userData->h_fft_result_im);
  free(userData->h_fft_result_re);
  free(userData->outputbuffer);
  free(userData->temp);
  free(userData->x_fft_result_im);
  free(userData->x_fft_result_re);
  free(userData->h);
  free(userData->conv);
  free(userData);

 cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  return 0;
}
