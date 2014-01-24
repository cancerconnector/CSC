// version for sorting out normal tissue consumption, will use small domain for computational efficiency 
// have to ensure that normal cells (0) consume at basal rate and then sort what the hell that is
// this version HAS NECROTIC CELLS, HAS NORMAL CELLS

import java.io.*;
import java.util.*;
import static java.lang.Math.*;

public class CA3 {
	public static String version="22nd Jan 2014";
    public boolean simulationFinished=false;
	MersenneTwisterFast random;
	public int mag = 1;
	public float [][] Oxygen=null; 
	public int [][] Age;
	int [][] Cells=null;  // Cells  and the cell type
						 // 0 -> No Cell
						 // 1 -> Stem Cells
						 // 2 -> Progenitor Cells
						 // 3 -> Mature/Differentiated cells
	int [][] Vasculature=null; // Initial vasculature as a boolean-like lattice
	Bag cellList=null; // Used in the iterateCells method. Defined here for performance issues
	boolean finishedRun=false;
	
	int size = 40; // Size of the system lattice
        int centre = size/2; //find centre of lattice
	int timestep=0; // Current Number of timesteps in the simulation
        int o2perTS = 1440; // 144000 Times we iterate oxygen per cell time step
        int dataWriteStep = 50; // how often we write down the data
        int dataReportStep = 50; // how often we write down the data

	// Model Parameters
	float initOxygen=0.0f;
	float[] proliferation = {
	 0.005f,  /* healthy cells */ 
	 0.005f, /* stem cells*/
	 0.005f, /* progenitor cells */
	 0.0f, /*mature/differentiated cells */
	 0.0f, /*necrotic area */ // added this cell type
	};
	
        float kDe = 0.1f; // 0.001728f;
        float consumptionBasal = 0.0001f; //********************** 2e-4 gives about diffusion length of 8, 1e-4 gives about 10
	float consumptionDivision = 4 * consumptionBasal; // The Oxygen consumption for cells that proliferate **was 5, changing to 4 on 18 december 2013,JGS
        float hypoxia = 0.1f; // Hypoxia threshold

	
	int maxProDivisions;
	int maxMatureCellAge;
	float asymmetricRatio;
	float pMotility=0.00f; // COMPLETELY ARBITRARY INDEED, REVISIT
	float densityVasculature=0.004f; // 1/250
	float avgVesselRad=2;
    boolean radiotherapy=false;
	
	// Measurement and STATISTICS
	int births=0;
	int deaths=0;
    int[][] stemBirthCounter=null;

	// Probability of stem cells producing progenitor cells given by oxygen concentration.
	// Probability of progenitor cells of producing differentiated cells: same

	public CA3 (float SCSymmetricDivBalance, int maxDivisions, float densityV)
	{
		int time = (int) System.currentTimeMillis();
		random = new MersenneTwisterFast (time);
		asymmetricRatio=SCSymmetricDivBalance;
		maxProDivisions=maxDivisions;
		maxMatureCellAge=1;
		densityVasculature=densityV;
		reset();
		resetVasculature();
	}
	
	public void reset ()
	{
		Oxygen = new float [size][size];
		Age = new int [size][size];
		stemBirthCounter = new int [size][size];
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				Oxygen[i][j]=initOxygen;
				Age[i][j]=0;
				stemBirthCounter[i][j]=0;
			}

		resetCells();
	}

    // *******initialize******* domain with vessels
	void resetVasculature ()
	{
		Vasculature = new int [size][size];	
		
		//for (int i=0;i<size;i++)
		    //for (int j=0;j<size;j++)
					    // 			if ((random.nextFloat()<=densityVasculature) && (Vasculature[i][j]==0)) Vasculature[i][j]=1;
			 Vasculature[centre][centre]=1;
	}
	
    // *******initialize******* domain with cells
	void resetCells ()
	{
		int radius=50;
		Cells = new int [size][size];
		
		// Let's set the empty space first.
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++) {
				Cells[i][j]=0;
			}

		Cells[centre][centre]=1; // initialize with cancer if you like
	}
    // to measure euclidean distance (not used)
	final int distance (int x, int y, int i, int j)
	{
		double dis=Math.sqrt((x-i)*(x-i)+(y-j)*(y-j));

		return (int)Math.round(dis);
	}

    // This method makes the boundaries no flux (not used yet) ************

	final int[] convertCoordinatesNoFlux(int x, int y)
       {
               // This method makes the boundaries no flux

               if (x < 0) x = 1;
               else if (x > size - 1) x = (size - 2);
               if (y < 0) y = 1;
               else if (y > size - 1) y = (size - 2);
               int[] result = new int[2];
               result[0] = x; result[1] = y;
               return result;
       }

     // This method makes the boundaries periodic

	final int[] convertCoordinates(int x, int y)
       {

               if (x < 0) x = size - 1;
               else if (x > size - 1) x = 0;
               if (y < 0) y = size - 1;
               else if (y > size - 1) y = 0;
               int[] result = new int[2];
               result[0] = x; result[1] = y;
               return result;
       }

	
	public void nextTimeStep ()
	{
	    // if (timestep==1000)   //zeros out vasculature at t=1000
            // for (int i=0;i<size;i++)
		// for (int j=0;j<size;j++) Vasculature[i][j]=0;
		births=0;
		deaths=0;
		for (int i=0;i<o2perTS;i++) iterateOxygen();
		iterateCells();
		radiotherapy=false; // Just for paranoia, remove!
		
		//NEW
		int totalCells=0;
		int totalHealthy=0;
		int totalStem=0;
		int totalProgenitor=0;
		int totalMature=0;
		for (int i=0;i<size;i++)
			for (int j=0;j<size;j++)
				if (Cells[i][j]<4) {
					totalCells++;
					if (Cells[i][j]==0) totalHealthy++;
					else if (Cells[i][j]==1) totalStem++;
					else if (Cells[i][j]==2) totalProgenitor++;
					else if (Cells[i][j]==3) totalMature++;
					else System.err.println ("wrong cell type");
				}
			
		// Not so new	
		if (timestep==0) System.out.println ("% Timestep\t Cells\t Stem Cells \t Progenitor\t Mature");
		if (timestep%dataReportStep==0) {		
		    System.out.println(timestep+"\t"+totalCells+"\t"+totalStem+"\t"+totalProgenitor+"\t"+totalMature+"\t"+births+"\t"+deaths+"\t"+((float)births/deaths));
		    System.err.println(asymmetricRatio+" "+maxProDivisions+" "+totalCells+" "+totalStem);}
		timestep++;
        
        // Finally let's write down the data
        if (timestep%dataWriteStep==0) {
            try {
                File dir = new File ("./text");
                dir.mkdir ();
		FileWriter outFile1 = new FileWriter("./text/cells"+timestep);
                PrintWriter outCells = new PrintWriter(outFile1);
                FileWriter outFile2 = new FileWriter("./text/oxygen"+timestep);
                PrintWriter outO2 = new PrintWriter(outFile2);
                // FileWriter outFile3 = new FileWriter("./text/stemBirths"+timestep);
                // PrintWriter outSB = new PrintWriter(outFile3);
                for (int i=0;i<size;i++) {
                    for (int j=0;j<size;j++){
			outCells.print(Cells[i][j]+", ");
                        outO2.print(Oxygen[i][j]+", ");
                        // outSB.print(stemBirthCounter[i][j]+", ");
                    }
                    outCells.println("");
                    outO2.println("");
                    // outSB.println("");
                }
		outCells.close();
                outO2.close();
                // outSB.close();
            } catch (IOException e) {
                e.printStackTrace();
                System.exit(-1);
            }
        }
	}
	

    // main CELL CA loop************
	public boolean iterateCells()
 	{
		if (cellList==null) cellList = new Bag (size*size); 
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++) {
			    if (Cells[i][j] < 4 )  { // All tumour cell types have Cell > 0, now 0 corresponds to 'healthy cells' that consume at basal rate only 
					int[] p = new int[2];
					p[0] = i; p[1] = j;
					cellList.add(p);
				}
			}

		while (cellList.size() != 0) {
			// Select the next lattice element at random
			int randomElemIndex=0;
			if (cellList.size()>1) randomElemIndex = random.nextInt(cellList.size()-1);
	   		int[] point = (int[])cellList.get(randomElemIndex); 
		   	int rI = point[0];
		   	int rJ = point[1];

			cellList.remove(randomElemIndex); // Remove it from the cell list
		   	int cell = Cells[rI][rJ];
	   
			// Cell death
			if ((Oxygen[rI][rJ]<hypoxia)) {
				if (random.nextInt(1)==0) { // 100% chances of dying in hypoxia. A BIT RANDOM!
					Age[rI][rJ]=0;
					Cells[rI][rJ]=4; // was 0, now making necrotic area (truly empty)
					deaths++;
					stemBirthCounter[rI][rJ]=0;
				}
			} else if ((cell==3) && (Age[rI][rJ]>100*maxMatureCellAge)) { // added 10* to allow for an pdate each celltimestep/10 **************************
				Age[rI][rJ]=0;
				Cells[rI][rJ]=4; // was 0, now making necrotic area (truly empty)
				deaths++;

			} else if ( (radiotherapy) && (cell==2) && (random.nextFloat()>Oxygen[rI][rJ]) ) {
                // Radiotherapy
                Age[rI][rJ]=0;
				Cells[rI][rJ]=4; // was 0, now making necrotic area (truly empty)
				deaths++;
            }
			// Do we have space for division or motility ?
			// else if (vacantSites(rI,rJ)>0) // this worked

			//healthy division
			else if ((cell==0) && (vacantSites(rI,rJ)>0)) {
			   if (proliferation[cell]>=random.nextFloat()) {// If tossing the coin we are to proliferate...
			       if (Oxygen[rI][rJ]>consumptionDivision) { // AND the oxygen concentration is enough for division...
				      Oxygen[rI][rJ]=Oxygen[rI][rJ]-consumptionDivision;
				      int[] daughter = findEmptySite (rI,rJ);
				births++;
				Cells[daughter[0]][daughter[1]]=0;
				}
			   }
		}
			// cancer division
			   else if ((vacantSitesCancer(rI,rJ)>0) && (cell>0))
			   if (proliferation[cell]>=random.nextFloat()) {// If tossing the coin we are to proliferate...
				if ((cell==1) || ((cell==2) && (Age[rI][rJ]<maxProDivisions))) // AND the cell is healthy, stem or progenitor ...
				    if (Oxygen[rI][rJ]>consumptionDivision) { // AND the oxygen concentration is enough for division...
							Oxygen[rI][rJ]=Oxygen[rI][rJ]-consumptionDivision;
                            int[] daughter = findEmptySiteCancer (rI,rJ); // changed to FESC from FES
                        if ((daughter[0]==0) || (daughter[0]==size) || (daughter[1]==0)||(daughter[1]==size)) {simulationFinished=true;}
							// to here
							births++; 
							if (cell==1) { // stem cell
							    stemBirthCounter[rI][rJ]++;
								if (Oxygen[rI][rJ]>hypoxia) {
                                    if (asymmetricRatio>random.nextFloat()) {Cells[daughter[0]][daughter[1]]=1;stemBirthCounter[daughter[0]][daughter[1]]=stemBirthCounter[rI][rJ];}
                                    else {Cells[daughter[0]][daughter[1]]=2;stemBirthCounter[daughter[0]][daughter[1]]=0;} // Otherwise differentiate
                                } else { // Only if there's hypoxia...
                                    float newASR=asymmetricRatio+0*(hypoxia-Oxygen[rI][rJ]); 
				    // jake messing about here with the ^^ multiplier. 
                                    if (newASR>random.nextFloat()) {Cells[daughter[0]][daughter[1]]=1;stemBirthCounter[daughter[0]][daughter[1]]=stemBirthCounter[rI][rJ];}
                                    else {Cells[daughter[0]][daughter[1]]=2;stemBirthCounter[daughter[0]][daughter[1]]=0;} // Otherwise differentiate
                                }
							} 
							else if (cell==2) {
								if (Age[rI][rJ]<maxProDivisions-1) {
									Cells[daughter[0]][daughter[1]]=2;
									Age[rI][rJ]++;
									Age[daughter[0]][daughter[1]]=Age[rI][rJ];
								} else {
									Cells[daughter[0]][daughter[1]]=3;
									Cells[rI][rJ]=3;
									Age[rI][rJ]=0;
									Age[daughter[0]][daughter[1]]=Age[rI][rJ];
								}
							}
						}
				} else if (pMotility>random.nextFloat()) { // Maybe we can try migration?
						int[] daughter = findEmptySite (rI,rJ);
						Cells[daughter[0]][daughter[1]]=cell;
						Cells[rI][rJ]=0;
						Age[daughter[0]][daughter[1]]=Age[rI][rJ];
						Age[rI][rJ]=0;
						System.err.println ("moving "+rI+", "+rJ);
				}
			// Aging for mature cells
			if (cell==3) Age[rI][rJ]++;
	  	}
 	return true;
 }
    // this method finds vacant sites for HEALTHY cells
final int vacantSites (int x, int y) 

{
	// We assume that necrotic material counts as neighbour? NOTE: MOORE NEIGHBORHOOD
	int total=0;
	int[] p = new int [2];
		
	p=convertCoordinates (x+1,y-1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x+1,y);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x+1,y+1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x,y-1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x,y+1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x-1,y-1);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x-1,y);
	if (Cells[p[0]][p[1]]==4) total++;
	p=convertCoordinates (x-1,y+1);
	if (Cells[p[0]][p[1]]==4) total++;
	return total;
}	

int[] findEmptySite (int x, int y)
{
	LinkedList vacantSites = new LinkedList();
	int[] tp1 = new int[2];
	int[] tp2 = new int[2];
	int[] tp3 = new int[2];
	int[] tp4 = new int[2];
	int[] tp5 = new int[2];
	int[] tp6 = new int[2];
	int[] tp7 = new int[2];
	int[] tp8 = new int[2];
	
	tp1=convertCoordinates (x+1,y-1);
	if (Cells[tp1[0]][tp1[1]]==4) vacantSites.add(tp1);
	tp2=convertCoordinates (x+1,y);
	if (Cells[tp2[0]][tp2[1]]==4) vacantSites.add(tp2);
	tp3=convertCoordinates (x+1,y+1);
	if (Cells[tp3[0]][tp3[1]]==4) vacantSites.add(tp3);
	tp4=convertCoordinates (x,y-1);
	if (Cells[tp4[0]][tp4[1]]==4) vacantSites.add(tp4);
	tp5=convertCoordinates (x,y+1);
	if (Cells[tp5[0]][tp5[1]]==4) vacantSites.add(tp5);
	tp6=convertCoordinates (x-1,y-1);
	if (Cells[tp6[0]][tp6[1]]==4) vacantSites.add(tp6);
	tp7=convertCoordinates (x-1,y);
	if (Cells[tp7[0]][tp7[1]]==4) vacantSites.add(tp7);
	tp8=convertCoordinates (x-1,y+1);
	if (Cells[tp8[0]][tp8[1]]==4) vacantSites.add(tp8);
	
	// Now let's see where.
	if (vacantSites.size() > 0) { // Now choose a vacant one, otherwise return the original location
		// pick a vacant site and return it
		int vacantElemIndex = random.nextInt(vacantSites.size());
		int[] p = (int[])vacantSites.get(vacantElemIndex);
		return (int[])p;	
	} else {
		int[] p = new int[2];
		p[0] = x; p[1] = y; // Just return the original
		System.out.println ("wrong!:"+vacantSites (x,y)+" - "+vacantSites.size());
		return p;
	}

}

   

    // a method for finding vacant sites for CANCER CELLS
final int vacantSitesCancer (int x, int y)
{
	// We assume that necrotic material and healthy cells DO NOT count as filling neighbouring sites. NOTE: MOORE NEIGHBORHOOD
	int total=0;
	int[] p = new int [2];
		
	p=convertCoordinates (x+1,y-1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++; //making into a logical OR to include possibility of cells dividing into areas of necrosis cells=4
	p=convertCoordinates (x+1,y);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinates (x+1,y+1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinates (x,y-1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinates (x,y+1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinates (x-1,y-1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinates (x-1,y);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	p=convertCoordinates (x-1,y+1);
	if ((Cells[p[0]][p[1]]==0) || (Cells[p[0]][p[1]]==4)) total++;
	return total;
	}

int[] findEmptySiteCancer (int x, int y)
{
	LinkedList vacantSitesCancer = new LinkedList();
	int[] tp1 = new int[2];
	int[] tp2 = new int[2];
	int[] tp3 = new int[2];
	int[] tp4 = new int[2];
	int[] tp5 = new int[2];
	int[] tp6 = new int[2];
	int[] tp7 = new int[2];
	int[] tp8 = new int[2];
	
	tp1=convertCoordinates (x+1,y-1);
	if ((Cells[tp1[0]][tp1[1]]==0) || (Cells[tp1[0]][tp1[1]]==4)) vacantSitesCancer.add(tp1);
	tp2=convertCoordinates (x+1,y);
	if ((Cells[tp2[0]][tp2[1]]==0) || (Cells[tp2[0]][tp2[1]]==4)) vacantSitesCancer.add(tp2);
	tp3=convertCoordinates (x+1,y+1);
	if ((Cells[tp3[0]][tp3[1]]==0) || (Cells[tp3[0]][tp3[1]]==4)) vacantSitesCancer.add(tp3);
	tp4=convertCoordinates (x,y-1);
	if ((Cells[tp4[0]][tp4[1]]==0) || (Cells[tp4[0]][tp4[1]]==4)) vacantSitesCancer.add(tp4);
	tp5=convertCoordinates (x,y+1);
	if ((Cells[tp5[0]][tp5[1]]==0) || (Cells[tp5[0]][tp5[1]]==4)) vacantSitesCancer.add(tp5);
	tp6=convertCoordinates (x-1,y-1);
	if ((Cells[tp6[0]][tp6[1]]==0) || (Cells[tp6[0]][tp6[1]]==4)) vacantSitesCancer.add(tp6);
	tp7=convertCoordinates (x-1,y);
	if ((Cells[tp7[0]][tp7[1]]==0) || (Cells[tp7[0]][tp7[1]]==4)) vacantSitesCancer.add(tp7);
	tp8=convertCoordinates (x-1,y+1);
	if ((Cells[tp8[0]][tp8[1]]==0) || (Cells[tp8[0]][tp8[1]]==4)) vacantSitesCancer.add(tp8);
	
	// Now let's see where.
	if (vacantSitesCancer.size() > 0) { // Now choose a vacant one, otherwise return the original location
		// pick a vacant site and return it
		int vacantElemIndexCancer = random.nextInt(vacantSitesCancer.size());
		int[] p = (int[])vacantSitesCancer.get(vacantElemIndexCancer);
		return (int[])p;	
	} else {
		int[] p = new int[2];
		p[0] = x; p[1] = y; // Just return the original
		System.out.println ("wrong!:"+vacantSitesCancer (x,y)+" - "+vacantSitesCancer.size());
		return p;
	}

}

    //////////////******************

    public void iterateOxygen()
    {
        //
        float[][] newOxygen = new float[size][size];
        for (int rI = 0; rI < size; rI++)
            for (int rJ = 0; rJ < size; rJ++) {
                // Determine the actual coordinates for top (-1,0), left(0,-1), right(0,1), below(1,0)
                // using periodic boundary conditions
                int[] top = convertCoordinatesNoFlux(rI - 1, rJ);
                int[] left = convertCoordinatesNoFlux(rI, rJ - 1);
                int[] right = convertCoordinatesNoFlux(rI, rJ + 1);
                int[] below = convertCoordinatesNoFlux(rI + 1, rJ);
                // Diffusion
                newOxygen[rI][rJ]
                    = Oxygen[rI][rJ] + (kDe *
                                    (Oxygen[top[0]][top[1]]
                 + Oxygen[left[0]][left[1]]
                 + Oxygen[right[0]][right[1]]
                 + Oxygen[below[0]][below[1]]
                 - 4.0f * Oxygen[rI][rJ]));

                // Consumption
                if (Cells[rI][rJ]<4) { // just changing this rule to <4 from >0 should allow all 0 positions to consume oxygen at the basal rate
		    newOxygen[rI][rJ] = newOxygen[rI][rJ] - consumptionBasal; ///o2perTS;  I have now rescaled consumptionBasal to be what the consumption should be per OXtimestep
                }
			
                // Production 
                if ((Vasculature[rI][rJ] == 1)) { //&& ((Cells[rI][rJ]==0) || (Cells[rI][rJ]==4))) {  // here is where the vascular occlusion takes place. 
                    newOxygen[rI][rJ] = 1.0f;
		}
                // Sanity check
		//	if (newOxygen[rI][rJ]>1.0f) {
		//   newOxygen[rI][rJ]=1.0f;
		//    System.out.println ("wrong!: Oxygen greater than one");}
		// else 
		    if (newOxygen[rI][rJ]<0.0f) {
			newOxygen[rI][rJ]=0.0f;}
		//System.out.println ("wrong!: Oxygen less than one");}
            }
        Oxygen = newOxygen;
    }
    

	public float [][] getOxygen() { return Oxygen; }
	public int [][] getCells () { return Cells; }
	public int [][] getVasculature() {return Vasculature;}
  
    ///////  this is to run CA without Vis ***********
    // public static void main(String args[]) {
     //        int maxTS=1000;
	//        CA3 ca;
	//		System.err.println ("# CA version:"+CA3.version);
		//        float SCSymmetricDivBalance=0.2f;
	//		int maxDivisions=10;
		//		float densityV=0.04f;
		//		if (args.length==4) {
		    //			SCSymmetricDivBalance = Float.parseFloat (args[0]);
			//			maxDivisions = Integer.parseInt (args[1]);
			//			maxTS=Integer.parseInt (args[2]);
			//			densityV=Float.parseFloat(args[3]);
			//			System.err.println ("Balance: "+SCSymmetricDivBalance+" maxDiv: "+maxDivisions+" maxTS: "+ maxTS);
			//		} else {
		    //			System.err.println ("Arguments needed: s/a maxDivisions timesteps, densityV");
			//			System.exit(-1);
			//		}
		//        ca = new CA3(SCSymmetricDivBalance, maxDivisions, densityV);
	//        for (int ts=0;ts<maxTS;ts++) ca.nextTimeStep();
	//    }
    ///////  this is to run CA without Vis ***********
};
