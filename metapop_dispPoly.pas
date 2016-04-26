//------------------------------------------------------------------------------
// metapop_dispPoly
//
// Emanuel A. Fronhofer
// november 2010
//
//------------------------------------------------------------------------------

// Copyright (C) 2015  Emanuel A. Fronhofer

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

program metapop_dispPoly;

uses Math;

//------------------------------------------------------------------------------
//constants --------------------------------------------------------------------

const
XMAX = 50;       	// no. of patches in the world (in x and y dim.)
YMAX = 20;

KMAX = 10000;       // habitat capacity
NO_ALLELES = 2;   	// no. of alleles per locus

RS = 1;           	// random seed
VARI = 0.1;       	// variance used for mutation

MAXRUNS = 25;		// no of repeats

MUT_DR = 0.0001;	// mutation rate for dominant/recessive


//------------------------------------------------------------------------------
//types ------------------------------------------------------------------------

type

 // One Individual ------------------------------------------------------------
  TIndividual = object
      locus_dispRate 	: array[1..NO_ALLELES] of real;                         // locus for dispersal rate
      locus_dispRate_DR : array[1.. NO_ALLELES] of boolean;						// is this allele dominant or recessive
      
      locus_tradeOff 	: array[1..NO_ALLELES] of real;                         // locus encoding reduction of lambda and my; relative to lambda (i.e. 0...1)
      locus_tradeOff_DR : array[1..NO_ALLELES] of boolean;                      // dominant/ recessive
  end;

  // One Patch -----------------------------------------------------------------
  TPatch = object
      males    : array[1..KMAX] of TIndividual;                                 // all males
      no_males : integer;                                                       // counts all males
      newMales : array[1..KMAX] of TIndividual;                                 // all new males (newly dispersed or born)
      no_newMales : integer;                                                    // counts all new males

      females  : array[1..KMAX] of TIndividual;                                 // all females
      no_females : integer;                                                     // counts all females
      newFemales : array[1..KMAX] of TIndividual;                               // all new females (newly dispersed or born)
      no_newFemales : integer;                                                  // counts all new females

  end;

  // World ---------------------------------------------------------------------
  TWorld = array[1..XMAX,1..YMAX] of TPatch;                                    //whole metapopulation


//------------------------------------------------------------------------------
// variables -------------------------------------------------------------------

var
  metapop : TWorld;             // the whole metapopulation

  N_null : integer;             // subpopulation size at simulation start
  occupancy_null : real;        // patch occupancy at simulation start
  capacity : integer;           // habitat capacity

  lambda_null : real;           // basic fertility i.e. mean no. of offspring
  beta : real;                  // konkurrenzkoeffizient
  alpha : real;					// allee effect
  sigma : real;                 // environmental fluctuations
  my_null : real;               // basic mortality rate
  gamma : integer;

  mutationRate : real;          // mutation rate
  p_extinction : real;          // extrinction probability
  
  generations : integer;        // no of generations
  time : integer;               // counter for generations

  no_emigrants : longint;       // no of emigrants

  data, input : text;           // data output and input files
  
  TradeOffCheckBox, 
  LinkageCheckBox, 
  TradeOffFunctionCheckBox, 
  NNDCheckBox, 
  MatingModeCheckBox : string;	// switch features on and off

  run_no : integer;				// counter for repeats


//------------------------------------------------------------------------------
// technical procedures --------------------------------------------------------

   // phenotype ----------------------------------------------------------------
   function phenotype(const allele_1, allele_2 : real; const allele_12, allele_22 : boolean):real;
   // calculates phenotype for both alleles (here: mean)
   // called by procedure "dispersal" and "reproduction"
   begin
   
   if (allele_12 = allele_22) then phenotype := (allele_1 + allele_2 ) / 2
   else
	begin
	//(TRUE = dominant)
	if (allele_12 = TRUE) then phenotype := allele_1;
	if (allele_22 = TRUE) then phenotype := allele_2;
	end;
   
   end;
   
   // analyse data ----------------------------------------------------------------
   procedure Analyse(const t:integer);
   var
     patchNo_x, patchNo_y : integer;                 	// counter for patch address
     f,m : integer;                              		// counter for females, males and alleles

   begin //1
     for patchNo_x :=1 to XMAX do
     for patchNo_y :=1 to YMAX do
      begin  //2
       // textfile output
       if ((t = 0) and (patchNo_x =1) and (patchNo_y = 1) and (run_no = 1)) then
          begin
          writeln(data,'      t','     sex','   dispRate_1','    d/r','   TradeOff_1','    d/r','   dispRate_2','    d/r','   TradeOff_2','    d/r', '     run');
          end;
       if (t > (0.9 * generations)) then
          begin          
          for f := 1 to metapop[patchNo_x,patchNo_y].no_females do writeln(data, t:7, '     f   ', metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[1]:7:5 ,'   ', metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[1],'   ', metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1]:7:5,'   ', metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1],'   ', metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[2]:7:5 ,'   ', metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[2],'   ', metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2]:7:5,'   ', metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2],'      ', run_no:5);          
          for m := 1 to metapop[patchNo_x,patchNo_y].no_males   do writeln(data, t:7, '     m   ', metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[1]:7:5   ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[1]  ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[1]:7:5  ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1]  ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[2]:7:5   ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[2]  ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[2]:7:5  ,'   ', metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2],'      ', run_no:5);               
          end;
       
       // end of output code ----------------------------------------------------
       end;   //2
    end;    //1


// functions -------------------------------------------------------------------
   // poisson distribution -----------------------------------------------------
   function Poisson(lamb:double):integer;
   begin
	// this function creates Poisson distributed random numbers
	// with mean lamb. the source code can be found in:
	// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
	// Numerical Recipes in Pascal. The Art of Scientific Computing. 
	// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
   end;

   // log Gauss distribution ---------------------------------------------------
   function loggauss(Fertility: real): real;               //geändert: sigma global!
   begin
	// calculates the mean random fertility
	// of a single female individual from a
	// Gaussian probability distribution,
	// using the variable "Fertility" as
	// mean value and "sigma" (environmental
	// fluctuations) as standard deviation. the source code can be found in:
	// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
	// Numerical Recipes in Pascal. The Art of Scientific Computing. 
	// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
   end;

   // mutation -----------------------------------------------------------------
   function mutate(const allele : real) : real;
   // introduces a mutation uniformely distributet around the mean "allele" with
   // variance VARI (see constants)!
   // from A. Poethke
   // called by procedure "Reproduction"
   var
      m1 : real;          
   begin
   m1 := allele + VARI * (2*random - 1);

   // values only between 0 and 1
   if m1 < 0 then m1 := abs(m1);
   if m1 > 1 then m1 := 2 - m1;

   mutate := m1;

   end;

   // incorporates trade off ---------------------------------------------------
   function lambdaTradeOff(const L, LR : real): real;
   // incorporates the trade off, i.e. reduces lambda according to trade off locus
   // linear function
   // called by function "lambda"
   // L : lambda null after incorporation of env. fluct.
   // LR : relative reduction of lambda
   begin
   lambdaTradeOff := L * (1 - LR)
   end;

   // lambda -------------------------------------------------------------------
   function lambda(const Nt, lambdaEF : real): real;
   // calculates the actual lambda (i.e. mean no of offspring) taking into account
   // logistic growth
   // called by procedure "reproduction"
   // Nt : actual pop size
   // lambdaRed is relative reduction of lambda
   var
      a, b : real;           	//susceptibility to crowding
      Nt1 : real;         		//pop  size in next generation

   begin
   // calculation susceptibility to crowding
   a := (lambda_null -1) / capacity; 
   
   // inroducing an allee effect
   b := power((Nt / capacity),2) / ( power((Nt / capacity),2) + power(alpha,2) );

   // calculate pop size in next generation with log. growth
   Nt1 := ( Nt * lambdaEF ) * b / power((1 + a * Nt), beta);

   if Nt > 0 then lambda := Nt1 / Nt
   else lambda := 0;

   end;

   // my -----------------------------------------------------------------------
   function my(const lambdaRed : real) : real;
   // calculated actual my (disp. mort.) after reduction by trade off
   // lambdaRed is relative reduction of lambda
   //called by procedure "dispersal"
   var
      my1 : real; 
   begin
  if TradeOffFunctionCheckBox = 'linear' then
   begin
   // linear function
   my1 := my_null * (1 - gamma * lambdaRed);

   // not negative
   if my1 < 0 then my := 0
   else my := my1;
   end;
	
   if TradeOffFunctionCheckBox = 'negexp' then
   begin
   // neg. exp. function
   // here "slopeMyCurve" is called gamma
   my := my_null * exp(-gamma * lambdaRed);
	end;
	
   end;

// biol. procedures ------------------------------------------------------------
   // new patch choice ---------------------------------------------------------
   procedure NewPatchChoice (const x,y : integer; var xn, yn : integer);
   var
   nnd_no : integer;    // target patch no for nnd
   begin   //1

   if NNDCheckBox = 'no' then
      begin     //2
        // new patch address if global dispersal
        repeat
              xn := random(XMAX) + 1;
              yn := random(YMAX) + 1;
              until (xn <> x) or (yn <> y);
      end     //2
   else
       begin    //2
       // new patch address if nnd
       nnd_no := random(8) + 1;
       case nnd_no of   //3
       1 : begin
           xn := x - 1;
           yn := y + 1;
           end;
       2:  begin
           xn := x;
           yn := yn +1;
           end;
       3:  begin
           xn := x +1;
           yn := y +1;
           end;
       4:  begin
           xn := x -1;
           yn := y;
           end;
       5:  begin
           xn := x +1;
           yn := y;
           end;
       6:  begin
           xn := x -1;
           yn := y-1;
           end;
       7:  begin
           xn := x;
           yn := y-1;
           end;
       8:  begin
           xn := x+1;
           yn := y-1;
           end;
       end;   //3

   if xn < 1 then xn := XMAX;
   if xn > XMAX then xn := 1;
   if yn < 1 then yn := YMAX;
   if yn > YMAX then yn := 1;

   end;   //2

   end;  //1

   // dispersal ----------------------------------------------------------------
   procedure dispersal;
   var
      patchNo_x, patchNo_y : integer;            // coordinates of actual patch
      f,m : integer;                             // counter for female and males
      newPatchNo_x, newPatchNo_y : integer;      // coordinates of patch to which an individual disperses
      help : integer;                            // counter

   begin    //1

   //reset counter for no of emigrants
   no_emigrants := 0;

   for patchNo_x := 1 to XMAX do                 //loops for every patch
   for patchNo_y := 1 to YMAX do
       begin   //2

       // if there are females...
       if metapop[patchNo_x,patchNo_y].no_females > 0 then
          begin //3

          f := 0;  // counter for no of females
          while f < metapop[patchNo_x,patchNo_y].no_females do      // loop for every female
              begin //4
              inc(f);
              // should they disperse: calculated with function "phenotype" (mean of both allels)
              if random < phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[1],metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[2], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[1],metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[2])  then
                 begin //5

                 // counter for no of emigrants
                 inc(no_emigrants);

                 //new patch is only reached if indiv survives!
                 // here: calculation of my (see function for reduction through trade-off)
                 if random > my(phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1],metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1],metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2])) then
                 begin  //6

                        // new patch address
                        NewPatchChoice(patchNo_x,patchNo_y,newPatchNo_x,newPatchNo_y);

                        // copy indiv. to new patch
                        inc(metapop[newPatchNo_x,newPatchNo_y].no_newFemales);
                        metapop[newPatchNo_x,newPatchNo_y].newFemales[metapop[newPatchNo_x,newPatchNo_y].no_newFemales] := metapop[patchNo_x,patchNo_y].females[f] ;
                 end;  //6

                 // erase idiv. from old patch
                 metapop[patchNo_x,patchNo_y].females[f] := metapop[patchNo_x,patchNo_y].females[metapop[patchNo_x,patchNo_y].no_females];
                 dec(metapop[patchNo_x,patchNo_y].no_females);
				dec(f);
                 end; //5
              end; // 4
          end;  //3

          // if there are males...
          if metapop[patchNo_x,patchNo_y].no_males > 0 then
          begin //3

          m := 0; //counter for males
          while m < metapop[patchNo_x,patchNo_y].no_males do      // loop for every male
              begin //4
              inc(m);

              // should they disperse: calculated with function "phenotype" (mean of both allels)
              if random < phenotype(metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[1],metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[2], metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[1],metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[2])  then
                 begin //5

                 // counter for no of emigrants
                 inc(no_emigrants);

                 //new patch is only reached if indiv survives!
                 // here: calculation of my (see function for reduction through trade-off)
                 if random > my(phenotype(metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[1],metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1],metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2])) then
                 begin  //6

                        // new patch address
                         NewPatchChoice(patchNo_x,patchNo_y,newPatchNo_x,newPatchNo_y);

                        // copy indiv. to new patch
                        inc(metapop[newPatchNo_x,newPatchNo_y].no_newMales);
                        metapop[newPatchNo_x,newPatchNo_y].newMales[metapop[newPatchNo_x,newPatchNo_y].no_newMales] := metapop[patchNo_x,patchNo_y].males[m] ;
                 end;  //6

                 // erase idiv. from old patch
                 metapop[patchNo_x,patchNo_y].males[m] := metapop[patchNo_x,patchNo_y].males[metapop[patchNo_x,patchNo_y].no_males];
                 dec(metapop[patchNo_x,patchNo_y].no_males);
					dec(m);
                 end; //5
              end; // 4
          end;  //3

       end;   //2

   //now that the dispersal phase is over, the residents and the dispersers are merged in one verctor!
   for patchNo_x := 1 to XMAX do           //loops for every patch
   for patchNo_y := 1 to YMAX do
       begin   //2

       // if there are females...
       if metapop[patchNo_x,patchNo_y].no_newFemales > 0 then
          begin //3

          help := metapop[patchNo_x,patchNo_y].no_females; // counter for array newFemales
          for f:= 1 to metapop[patchNo_x,patchNo_y].no_newFemales do      // loop for every female
              begin //4
              inc(help);
              // insert new  females in resident female array
              metapop[patchNo_x,patchNo_y].females[help] := metapop[patchNo_x,patchNo_y].newFemales[f];
              end; //4

          // increase counter for no of females
          metapop[patchNo_x,patchNo_y].no_females := metapop[patchNo_x,patchNo_y].no_females + metapop[patchNo_x,patchNo_y].no_newFemales;

          //reset counter for new indiv
          metapop[patchNo_x,patchNo_y].no_newFemales := 0;

          end; //3

       // if there are males...
       if metapop[patchNo_x,patchNo_y].no_newMales > 0 then
          begin //3

          help := metapop[patchNo_x,patchNo_y].no_males; // counter for array newMales
          for m:= 1 to metapop[patchNo_x,patchNo_y].no_newMales do      // loop for every male
              begin //4
              inc(help);
              // insert new  males in resident male array
              metapop[patchNo_x,patchNo_y].males[help] := metapop[patchNo_x,patchNo_y].newMales[m];
              end; //4

          // increase counter for no of females
          metapop[patchNo_x,patchNo_y].no_males := metapop[patchNo_x,patchNo_y].no_males + metapop[patchNo_x,patchNo_y].no_newMales;

          //reset counter for new indiv
          metapop[patchNo_x,patchNo_y].no_newMales := 0;

          end; //3

       end; //2

   end; //1

   // reproduction -------------------------------------------------------------
   procedure reproduction;
    var
      patchNo_x, patchNo_y : integer;       // coordinates of actual patch
      f,m,o : integer;                      // counter for females, males and offspring
      no_offspring : integer;               // no. of offspring
      pop_size : longint;                   // actual pop. size
      allele_no : integer;                  // which allele is inherited
      lambdaEF : real;						// lambda
      delta, delta_old : real;				// distance in phenotype space between female and potential mate
      potmate : integer;					// no of potential mate when assortative mating is simulated

   begin //1

   for patchNo_x := 1 to XMAX do             //loops for every patch
   for patchNo_y := 1 to YMAX do
       begin   //2

       //reset counter for new indiv
       metapop[patchNo_x,patchNo_y].no_newFemales := 0;
       metapop[patchNo_x,patchNo_y].no_newMales := 0;

       //check whether there are females and males
       if (metapop[patchNo_x,patchNo_y].no_females > 0) and (metapop[patchNo_x,patchNo_y].no_males > 0) then
          begin //3

              // calculate actual pop. size
              pop_size := metapop[patchNo_x,patchNo_y].no_males + metapop[patchNo_x,patchNo_y].no_females;

          // calculate lambda after incorporation of var sigma (environmental fluctuation)
          if lambda_null > 0 then lambdaEF:= lambda(pop_size, loggauss(lambda_null))   //geändert: sigma global!
          else lambdaEF := 0;

          //loop for every female
          for f:=1 to metapop[patchNo_x,patchNo_y].no_females do
              begin //4

             //randomly choose a male for mating
              if MatingModeCheckBox = 'random' then m := random(metapop[patchNo_x,patchNo_y].no_males) + 1;
              
              // assortative mating
              if MatingModeCheckBox = 'assortative' then
              begin //5
              // initialize phenotypic distance (set to maximal)
              // start with randomly choosing a mate
              m := random(metapop[patchNo_x,patchNo_y].no_males) + 1;
              
              // calculate genetic distance to this randomly choosen mate
              delta_old := sqrt(  power( (phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2]) - phenotype(metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[1], metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1], metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2]) ),2)
							    + power( (phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[1], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[2], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[1], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[2]) - phenotype(metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[1], metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[2], metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[1], metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[2]) ),2)
							   );
                          
              for potmate := 1 to metapop[patchNo_x,patchNo_y].no_males do
				begin //6
				// calculate new genetic distance between the focal female and a potential new mate in phenotype space				
                  delta := sqrt(  power( (phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2]) - phenotype(metapop[patchNo_x,patchNo_y].males[potmate].locus_tradeOff[1], metapop[patchNo_x,patchNo_y].males[potmate].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].males[potmate].locus_tradeOff_DR[1], metapop[patchNo_x,patchNo_y].males[potmate].locus_tradeOff_DR[2]) ),2)
							    + power( (phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[1], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[2], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[1], metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[2]) - phenotype(metapop[patchNo_x,patchNo_y].males[potmate].locus_dispRate[1], metapop[patchNo_x,patchNo_y].males[potmate].locus_dispRate[2], metapop[patchNo_x,patchNo_y].males[potmate].locus_dispRate_DR[1], metapop[patchNo_x,patchNo_y].males[potmate].locus_dispRate_DR[2]) ),2)
							   );
				
				if (delta < delta_old) then
				begin //7
					m:= potmate;
					delta_old := delta;
				end; //7
				
				end; //6
			   end; //5


              // calculate no. of offspring
              // here the trade-off alleles of the female are sent to function lambda in order to calculate lambda
              // this function incorporates logistic growth
              // then the no. of offspring is calculated from a poisson distribution
              no_offspring := poisson(2 * lambdaTradeOff(lambdaEF, phenotype(metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1],metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2], metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1],metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2])) );
              
              if no_offspring > 0 then
                 begin //5
                 for o := 1 to no_offspring do
                     begin //6

                     //randomly assign offspring sex
                     if random < 0.5 then
                        //female offspring
                        begin // 7
                        inc(metapop[patchNo_x,patchNo_y].no_newFemales);

                        // determine genetics of offspring
                        // mutation rate is introduced
                        // mutation at trade off locus only if trade off is on!

                        //father's allele
                        //dispersal rate locus
                        // the allele is chosen randomly
                        allele_no := random(2)+1;       // either 1 or 2

                        metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate[1] := metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[allele_no];
                        if random < mutationRate then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate[1] := mutate(metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate[1]);

						// dispersal rate locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[1] := metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[allele_no];
						if random < MUT_DR then 
							begin
							if metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[1] = TRUE then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[1] := FALSE
							else metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[1] := TRUE;
							end;

                        // trade off locus
                        // either both loci are linked, i.e. allele_no is the same, or not
                        if LinkageCheckBox = 'no' then allele_no := random(2)+1;

                        metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff[1] := metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[allele_no];
                        if (random < mutationRate) and (TradeOffCheckBox = 'yes') then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff[1] := mutate(metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff[1]);

						// trade off locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[1] := metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[allele_no];
						if random < MUT_DR then
						begin
						if metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[1] = TRUE then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[1] := FALSE
						else metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[1] := TRUE;
						end;

                        // mother's allele
                        // dispersal rate locus
                        // the allele is chosen randomly
                        allele_no := random(2)+1;       // either 1 or 2

                        metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate[2] := metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[allele_no];
                        if random < mutationRate then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate[2] := mutate(metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate[2]);

						// dispersal rate locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[2] := metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[allele_no];
						if random < MUT_DR then
						begin
						if metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[2] = TRUE then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[2] := FALSE
						else metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_dispRate_DR[2] := TRUE;
						end;

                        // trade off locus
                        // either both loci are linked, i.e. allele_no is the same, or not
                        if LinkageCheckBox = 'no' then allele_no := random(2)+1;

                        metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff[2] := metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[allele_no];
                        if (random < mutationRate) and (TradeOffCheckBox = 'yes') then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff[2] := mutate(metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff[2]);

						// trade off locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[2] := metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[allele_no];
						if random < MUT_DR then
						begin
						if metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[2] = TRUE then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[2] := FALSE
						else metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_tradeOff_DR[2] := TRUE;
						end;

                        end   // 7
                     else
                        // male offspring
                        begin // 7
                        inc(metapop[patchNo_x,patchNo_y].no_newMales);

                        //determine genetics of offspring
                        // mutation rate is introduced

                        //father's allele
                        //dispersal rate locus
                        // the allele is chosen randomly
                        allele_no := random(2)+1;       // either 1 or 2

                        metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate[1] := metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[allele_no];
                        if random < mutationRate then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate[1] := mutate(metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate[1]);

						// dispersal rate locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[1] := metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[allele_no];
						if random < MUT_DR then
						begin
						if metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[1] = TRUE then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[1] := FALSE
						else metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[1] := TRUE;
						end;

                        //trade off locus
                        // either both loci are linked, i.e. allele_no is the same, or not
                        if LinkageCheckBox = 'no' then allele_no := random(2)+1;

                        metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff[1] := metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[allele_no];
                        if (random < mutationRate) and (TradeOffCheckBox = 'yes') then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff[1] := mutate(metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff[1]);

						// trade off locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[1] := metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[allele_no];
						if random < MUT_DR then
						begin
						if metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[1] = TRUE then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[1] := FALSE
						else metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[1] := TRUE;
						end;

                        //mother's allele (chosen randomly)
                        //dispersal rate locus
                        // the allele is chosen randomly
                        allele_no := random(2)+1;       // either 1 or 2

                        metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate[2] := metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[allele_no];
                        if random < mutationRate then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate[2] := mutate(metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate[2]);

						// dispersal rate locus dominant or recessive with possible mutation
						metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[2] := metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[allele_no];
						if random < MUT_DR then
						begin
						if metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[2] = TRUE then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[2] := FALSE
						else metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_dispRate_DR[2] := TRUE;
						end;

                        // trade off locus
                        // either both loci are linked, i.e. allele_no is the same, or not
                        if LinkageCheckBox = 'no' then allele_no := random(2)+1;

                        metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff[2] := metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[allele_no];
                        if (random < mutationRate) and (TradeOffCheckBox = 'yes') then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff[2] := mutate(metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff[2]);
                        
                        // trade off locus dominant or recessive with possible mutation
                        metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[2] := metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[allele_no];
                        if random < MUT_DR then
                        begin
                        if metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[2] = TRUE then metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[2] := FALSE
                        else metapop[patchNo_x,patchNo_y].newMales[metapop[patchNo_x,patchNo_y].no_newMales].locus_tradeOff_DR[2] := TRUE;
                        end;
                        
                        
                        end;  // 7
                     end; //6
                 end;  //5
              end; //4
          end;  //3
       end; //2
   end; //1

   // death of annual organisms -------------------------------------------------
   procedure death;
   var
     patchNo_x, patchNo_y : integer;     // coordinates of actual patch
     m,f : integer;                      // counter for males and females

   begin // 1

   for patchNo_x := 1 to XMAX do         //loops for every patch
   for patchNo_y := 1 to YMAX do
       begin   //2

       // replace adult female vektor with offspring vector
       for f:= 1 to metapop[patchNo_x,patchNo_y].no_newFemales do
           begin //3
           metapop[patchNo_x,patchNo_y].females[f] := metapop[patchNo_x,patchNo_y].newFemales[f];
           end;  //3

       //reset female counters
       metapop[patchNo_x,patchNo_y].no_females := metapop[patchNo_x,patchNo_y].no_newFemales;
       metapop[patchNo_x,patchNo_y].no_newFemales := 0;

       // replace adult male vektor with offspring vector
       for m:= 1 to metapop[patchNo_x,patchNo_y].no_newMales do
           begin //3
           metapop[patchNo_x,patchNo_y].males[m] := metapop[patchNo_x,patchNo_y].newMales[m];
           end;  //3

       //reset male counters
       metapop[patchNo_x,patchNo_y].no_males := metapop[patchNo_x,patchNo_y].no_newMales;
       metapop[patchNo_x,patchNo_y].no_newMales := 0;

       // random patch extinction
       if random < p_extinction then
       begin //3
       metapop[patchNo_x,patchNo_y].no_females := 0;
       metapop[patchNo_x,patchNo_y].no_males := 0;
       end;  //3

       end; //2
   end; // 1

// simulation ---------------------------------------------------------------------
   // initialize --------------------------------------------------------
   procedure Initialize;
   var
      patchNo_x, patchNo_y : integer;
      f, m : integer;
             outputfilename : string;	

   begin //1      

        // parameter einlesen
        assign(input, 'input.in');
        reset(input);
        
        readln(input); readln(input); readln(input); readln(input);
        readln(input, occupancy_null); readln(input);
		readln(input, N_null); readln(input);       
		readln(input, capacity); readln(input);
		readln(input, beta); readln(input);
		readln(input, lambda_null); readln(input);
		readln(input, sigma); readln(input);
		readln(input, my_null); readln(input);
		readln(input, gamma); readln(input);
		readln(input, mutationRate); readln(input);
		readln(input, p_extinction); readln(input);
		readln(input, generations); readln(input); 
		readln(input, alpha); readln(input); readln(input);
		readln(input, TradeOffCheckBox); readln(input);
		readln(input, LinkageCheckBox); readln(input);
		readln(input, NNDCheckBox); readln(input);
		readln(input, TradeOffFunctionCheckBox); readln(input);
		readln(input, MatingModeCheckBox);
		
		close(input);
		
		if (run_no = 1) then
		begin
        outputfilename := 'results.out';
        assign(data, outputfilename);
        rewrite(data);
		end;
        
        //initialize patches in world
        for patchNo_x := 1 to XMAX do   // loop in x dir.
        for patchNo_y := 1 to YMAX do  // loop in y dir.
        begin   //2
             // initialize all patches empty
             metapop[patchNo_x,patchNo_y].no_females := 0;
             metapop[patchNo_x,patchNo_y].no_newFemales := 0;

             metapop[patchNo_x,patchNo_y].no_males := 0;
             metapop[patchNo_x,patchNo_y].no_newMales := 0;

             // takes occupancy(t0) into account
             if random < occupancy_null then
             begin //3
                   //initialize females
                   metapop[patchNo_x,patchNo_y].no_females := round(N_null / 2);

                   for f := 1 to metapop[patchNo_x,patchNo_y].no_females do
                   begin //4
                         // random initialization for dispersal rate locus
                         metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[1] := random;
                         metapop[patchNo_x,patchNo_y].females[f].locus_dispRate[2] := random;
                         
                         // random initialization of dominance/ recessivity (TRUE = dominant)
                         if random < 0.5 then metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[1] := FALSE
							else metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[1] := TRUE;
                         if random < 0.5 then metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[2] := FALSE
							else metapop[patchNo_x,patchNo_y].females[f].locus_dispRate_DR[2] := TRUE;

                         //random initialization for trade off locus if check box checked
                         if TradeOffCheckBox = 'yes' then
                         begin  //5
                              metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1] := random;
                              metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2] := random;
                              
                              // random initialization of dominance/ recessivity (TRUE = dominant)
                               if random < 0.5 then metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1] := FALSE
								  else metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1] := TRUE;
							   if random < 0.5 then metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2] := FALSE
								  else metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2] := TRUE;
                         end   //5
                         else
                         begin  //5
                              metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[1] := 0;
                              metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff[2] := 0;
                              
                              // random initialization of dominance/ recessivity (TRUE = dominant)
                              if random < 0.5 then metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1] := FALSE
								 else metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[1] := FALSE;
							  if random < 0.5 then metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2] := FALSE
								 else metapop[patchNo_x,patchNo_y].females[f].locus_tradeOff_DR[2] := FALSE;
                         end;  //5

                   end; //4

                   //initialize males
                   metapop[patchNo_x,patchNo_y].no_males := N_null - metapop[patchNo_x,patchNo_y].no_females;

                   for m := 1 to metapop[patchNo_x,patchNo_y].no_males do
                   begin //4
                         // random initialization for dispersal rate locus
                         metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[1] := random;
                         metapop[patchNo_x,patchNo_y].males[m].locus_dispRate[2] := random;
                         
                         // random initialization of dominance/ recessivity (TRUE = dominant)
                         if random < 0.5 then metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[1] := FALSE
							else metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[1] := TRUE;
                         if random < 0.5 then metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[2] := FALSE
							else metapop[patchNo_x,patchNo_y].males[m].locus_dispRate_DR[2] := TRUE;

                         // random initialization for trade off locus if check box checked
                         if TradeOffCheckBox = 'yes' then
                         begin //5
								metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[1] := random;
								metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[2] := random;
                         
							   // random initialization of dominance/ recessivity (TRUE = dominant)
							   if random < 0.5 then metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1] := FALSE
								  else metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1] := TRUE;
							   if random < 0.5 then metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2] := FALSE
								  else metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2] := TRUE;                         
                         end  //5
                         else
                         begin //5
							   metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[1] := 0;
							   metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff[2] := 0;
                         
                               // random initialization of dominance/ recessivity (TRUE = dominant)
                               if random < 0.5 then metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1] := FALSE
								 else metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[1] := FALSE;
							   if random < 0.5 then metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2] := FALSE
								 else metapop[patchNo_x,patchNo_y].males[m].locus_tradeOff_DR[2] := FALSE;
                         end;  //5
                   end; //4

             end; //3
        end;  //2

        //analyse after initialization (writes first line to output file)
        Analyse(0);

   end;  //1

//----------------------------------------------------------------------
// main ----------------------------------------------------------------   
begin

	// intialize random number generator
    randseed := RS;  

	for run_no:=1 to MAXRUNS do
	begin     
	Initialize;										// initializes world and reads form input file
   
	for time:= 1 to generations do           		// time loop
       begin
       dispersal;                          			// juvenile dispersal
       reproduction;                       			// reproduction in new patch
       death;                              			// death of all adults
       
       //analyse only when output is needed
       if time = generations then Analyse(time);    // analysis and textfile output
       end;
	
	end;											// end of time loop
	close(data);
end.

