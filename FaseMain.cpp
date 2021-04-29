#include "Common.h"
#include "Fase.h"
#include "DynamicGraph.h"
#include "GraphMatrix.h"
#include "GraphUtils.h"
#include "Timer.h"
#include "Graph.h"
#include "Random.h"
#include "Stream.h"
#include "RecalculateStream.h"
#include "SubgraphStream.h"
#include "DESUStream.h"

using namespace std;

Graph *G;
Stream *S;
int K = 0;
double sampProb[MAXMOTIF], prob;
bool dir = false, detailed = false, draw = false, samp = false, largeScale = false, stream = false, weighted = false;
char ifilename [200];
char ofilename [200];
char sfilename [200];
FILE *outFile, *streamFile;
time_t t_start, t_end;
int zeroBased = 1;
int batchSize = 1;
bool debugprints = false;

void init()
{
  printf("  88888888888           ad88888ba   88888888888  \n"
         "  88                   d8\"     \"8b  88           \n"
         "  88                   Y8,          88           \n"
         "  88aaaaa  ,adPPYYba,  `Y8aaaaa,    88aaaaa      \n"
         "  88\"\"\"\"\"  \"\"     `Y8    `\"\"\"\"\"8b,  88\"\"\"\"\"      \n"
         "  88       ,adPPPPP88          `8b  88           \n"
         "  88       88,    ,88  Y8a     a8P  88           \n"
         "  88       `\"8bbdP\"Y8   \"Y88888P\"   88888888888  \n\n"
         "\tVersion: 1.0\n"
         "FaSE - Fast Subgraph Enumeration (with Sampling)\n"
         "\n\n\tPedro {Paredes, Ribeiro} - DCC/FCUP\n\n\n\n\n");
  t_start = time(0);
}

void displayHelp()
{
  printf("------------ FaSE Usage --------------\nMain Settings: ./FASE -s <Subgraph Size> -i <input file> [arguments...]\n\n\tAll commands:\n-h : Displays this help information\n-s <Integer> : Subgraph Size\n-i <Filename> : Name of input file (Format in Readme.txt)\n-d : Directed Subgraph (Default undirected)\n-o : Name of output file (Default is stdout)\n-S : Stream mode from console input\n-Sf : Stream updates with file\n-dt : Detailed Result (Displays all subgraph types and occurrences)\n-ls : Use a large scale representation (default is adjacency matrix)\n-z : Use 0-based input (Suitable for input files starting at node 0)\n-p <P1> <P2> ... <Ps> : Sets the sampling probabilities by depth (note that -s must have been selected first)\n-q : Ignore arguments and prompt input\n--------------------------------------\n");
}

void read(int argc, char **argv)
{
  int E, V, i, check = 0, itera = 0;
  ofilename[0] = '0';
  ofilename[1] = '\0';
  sfilename[0] = '0';
  sfilename[1] = '\0';
  S = new SubgraphStream();
  for (i = 1; i < argc; i++)
  {
    if (argv[i][0] != '-')
      continue;
    if (argv[i][1] == 'h')
    {
      displayHelp();
      K = 0;
      return;
    }

    if (argv[i][1] == 'd' && argv[i][2] == 't')
      detailed = true;
    else if (argv[i][1] == 'd' && argv[i][2] == 'r')
      draw = true;
    else if (argv[i][1] == 'd')
      dir = true;

    if (argv[i][1] == 'l' && argv[i][2] == 's')
      largeScale = true;

    if (argv[i][1] == 'z')
      zeroBased = 0;

    if(argv[i][1] == 'x')
      debugprints = true;

    if (argv[i][1] == 'w')
      weighted = true;

    if (argv[i][1] == 'i')
    {
      i++;
      strcpy(ifilename, argv[i]);
      check |= (1 << 0);
      continue;
    }
    if (argv[i][1] == 'S' && argv[i][2] == 'f'){
      //printf("stream file");
      i++;
      stream = true;
      strcpy(sfilename,argv[i]);
      continue;
    }
    else if(argv[i][1] == 'S'){
      stream = true;
    }

    
    if (argv[i][1] == 'R' && argv[i][2] == 's'){
      stream = true;
      delete S;
      S = new RecalculateStream();
    }

    if(argv[i][1] == 'B' && argv[i][2] == 's'){
      i++;
      batchSize = argv[i][0] - '0';
      int j = 1;
      while (argv[i][j] != '\0')
      {
        batchSize *= 10;
        batchSize += argv[i][j] - '0';
        j++;
      }
      
      continue;
    }else if(argv[i][1] == 'B'){

      delete S;
      S = new DESUStream();
      continue;
    }

    

    if (argv[i][1] == 's')
    {
      i++;
      K = argv[i][0] - '0';
      int j = 1;
      while (argv[i][j] != '\0')
      {
      	K *= 10;
      	K += argv[i][j] - '0';
      	j++;
      }
      check |= (1 << 1);
      continue;
    }

    if (argv[i][1] == 'o')
    {
      i++;
      strcpy(ofilename, argv[i]);
      continue;
    }

    if (argv[i][1] == 'p')
    {
      int j;
      for (j = 0, i++; j < K; j++, i++)
        sampProb[j] = atof(argv[i]);
      samp = true;
      continue;
    }

    if (argv[i][1] == 'q')
    {
      itera = 1;
      break;
    }
  }

  if (largeScale)
    G = new DynamicGraph();
  else
    G = new GraphMatrix();

  if (!itera)
  {
    if (check != (1 << 2) - 1)
    {
      K = 0;
      if (check != 0)
	printf("Warning: Incorrect number of necessary arguments provided\n");
      displayHelp();
      return;
    }

    GraphUtils::readFileTxt(G, ifilename, dir, weighted, zeroBased);
    G->sortNeighbours();
    G->makeArrayNeighbours();
    //if(G->matrixNeighbours()) printf("I created it here...\n");
    
   
    
    return;
  }

  // Direction
  printf("Directed? (Y/n) ");
  char chdir;
  scanf(" %c", &chdir);
  if (chdir == 'n' || chdir == 'N')
    dir = false;

  // Initial
  printf("Input 0 or 1 based: ");
  scanf("%d", &zeroBased);

  // Input filename
  printf("Insert input file name: ");
  scanf(" %s", ifilename);
  GraphUtils::readFileTxt(G, ifilename, dir, false, zeroBased);
  G->sortNeighbours();
  G->makeArrayNeighbours();

  // Subgraph Size
  printf("Input the value K of the subgraph search: ");
  scanf("%d", &K);

  // Input filename
  printf("Insert output file name or 0 to stdout: ");
  scanf(" %s", ofilename);
  if (ofilename[0] == '0' && ofilename[1] == '\0')
    outFile = stdout;
  else
    outFile = fopen(ofilename, "w");

  // Default Sampling probabilities
  if (!samp)
    for (i = 0; i < K; i++)
      sampProb[i] = 1.0;
}

void initSamplingProbabilities(Fase* fase)
{
  int i;
  prob = 1.0;
  for (i = 0; i < K; i++)
    prob *= sampProb[i];

  if (samp && fabs(prob) > 10e-7)
    fase->initSampling(K, sampProb);
  else
    prob = 1.0;
}



int main(int argc, char **argv)
{
  init();
  //printf("oi\n");
  read(argc, argv);
  //printf("in main\n");
  if (K <= 2 || K >= MAXMOTIF)
  {
    fprintf(stderr, "Subgraph size needs to be between 3 and %d...\n", MAXMOTIF - 1);
    return 1;
  }
  /*
  Random::init(time(NULL));
  Fase* fase = new Fase(G, dir);
  initSamplingProbabilities(fase);

  Timer::start();
  fase->runCensus(K);
  Timer::stop();

  output(fase);
  finish(fase);
  */
  
  S->createStream(G,samp,sampProb,dir,sfilename,ofilename,stream,zeroBased);
  
  S->runCensus(K,batchSize,debugprints);
  delete S;
  return 0;
}
