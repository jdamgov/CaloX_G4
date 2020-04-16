#ifndef CaloDataStruc_h
#define CaloDataStruc_h 1

struct CaloStepData {
   int pid;
   int trackid;
   double x;
   double y;
   double z;
   double steplength;
   double globaltime;
   double edep;

};

#endif
