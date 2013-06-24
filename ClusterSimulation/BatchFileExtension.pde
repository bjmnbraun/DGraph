/**
 * A few useful tools for automating dGraph, notably:
 * 
 * 1) Running the simulation multiple times from 1 configuration
 * 2) Automatically rotating / flipping the outputs so they 
 * "match" the best
 * 3) Automatically saving the transformed coordinates in 
 * easy-to-analyze formatted files in a specified directory
 * 
 **/
void doCommandLine(){
  //parse arguments
  String help = System.getProperty("help");
  if (help!=null){
    doHelp();
    thenExit();
  }
  boolean doAutomation = false;
  String numRuns = System.getProperty("runs");
  String outputDir = System.getProperty("out");
  String filePref = System.getProperty("prefix");
  String numIters = System.getProperty("iters");
  String params = System.getProperty("params");
  String inFasta = System.getProperty("in");
  String inCustomMatrix = System.getProperty("inCustom");
  String doPdf = System.getProperty("pdf");
  String zoom = System.getProperty("zoom");
  doAutomation |= numRuns!=null;
  doAutomation |= outputDir!=null;
  doAutomation |= filePref!=null;
  doAutomation |= numIters!=null;
  doAutomation |= inFasta!=null;
  doAutomation |= inCustomMatrix!=null;
  if (doAutomation){
    //Necessary params:
    if (outputDir==null){
      batchError("-Dout must be specified.");
    }
    if (inFasta==null){
      batchError("-Din must be specified.");
    }
    //Necessary, if not doing PDF
    if (doPdf==null){ //Pdf render doesn't require all this stuff.
      if (numRuns==null){
        batchError("-Druns must be specified.");
      }
      if (params==null){
        batchError("-Dparams must be specified.");
      }
      if (numIters==null){
        batchError("-Diters must be specified.");
      }
    } else {
      numRuns = "0";
      numIters = "0";
      params = null; //Use defaults
    }
    //Has a default / Optional:
    if (filePref==null){
      filePref="dGbatch";
    }
    if (zoom!=null){
      VIEW_SCALE_USER=max(1,new Integer(zoom));
    }
    //Ok, do it.
    runScript(new Integer(numRuns),new Integer(numIters),params,outputDir,filePref,inFasta,inCustomMatrix,doPdf!=null);
  }
  //If we get here, we didn't input any command line arguments.
  System.out.println("Dgraph can be automated: try adding -Dhelp to the arguments when running the program.");
}
private void runScript(int numTimes, int numIters, String params1, String outputDirectory, String filePrefix, String inFasta, String customMatrix, boolean pdfOutput){
  iS = true; //inScript flag.

  File outputDir = new File(outputDirectory).getAbsoluteFile();
  if (!outputDir.isDirectory() || !outputDir.exists()){
    batchError("-Dout must be a pre-existing directory."+outputDir);
  }
  File inputFasta = new File(inFasta).getAbsoluteFile();
  if (!inputFasta.exists()){
    batchError("-Din referred to bad/no file\n"+inFasta);
  }
  toRead = inFasta; //Read by clusterSimulation.
  if (customMatrix!=null){
    File inputCustomMatrix = new File(customMatrix);
    if (!inputCustomMatrix.exists()){
      batchError("-DinCustom referred to bad/no file\n"+inputCustomMatrix);
    }
    GetCustomDistances(customMatrix);//Read by clusterSimulation.
  }

  if (params1!=null){ //Valid if just doing PDF
    String[] params = params1.split(",");
    CURRENT_CUTOFF = new Float(params[0]);
    SHOW_LINES = new Float(params[1]);
    REMOVAL_ISLAND_FINDER = new Float(params[2]);
    timeStep = new Float(params[3]);
    fluidFriction = new Float(params[4]);
    mass = new Float(params[5]);
  }

  //First time, reinit.
  reinit();

  if (pdfOutput){
    pdfUrl = outputDir.toString()+File.separator+filePrefix+".pdf";
    singleFramePdf = true;
    return;
  }

  System.out.printf("Setup completed. Running simulations on %d nodes:",nodes.length);
  System.out.println();
  String[] filesWrote = new String[numTimes];
  float[] BestRuns = new float[numTimes];
  for(int run = 0; run < numTimes; run++){
    for(int iter = 0; iter < numIters; iter++){
      runActualSimulation();
    }
    //Write the result:
    String outFasta = outputDir.toString()+File.separator+filePrefix+run;
    saveNodeStates(outFasta+".fasta");
    filesWrote[run]=outFasta; //For analysis later.
    BestRuns[run] = TOTAL_VARIANCE;
    //Notify.
    System.out.print("A");
    if (run%20==19){
      System.out.println();
    }
    //Latter times, only scatter.
    scatterNodes();
  }
  System.out.println();
  System.out.println("Analyzing...");
  try {
    PrintWriter ana = new PrintWriter(new FileWriter(new File(outputDir.toString()+File.separator+filePrefix+".dginfo")));
    float bestGood = Float.MAX_VALUE; 
    int bestRunIndex = 0;
    ana.println("Start of figure errors:");
    for(int run = 0; run < numTimes; run++){
      float goodness = BestRuns[run];
      ana.printf("%10d%10.3f",run,goodness);
      ana.println();
      if (goodness < bestGood){
        bestGood = goodness;
        bestRunIndex = run;
      }
    }
    ana.println();
    ana.println("Best run: "+bestRunIndex);
    int numNodes = nodes.length; //Still a leftover.
    float[][] Solution1 = new float[numNodes][2];
    getCoordsFile(Solution1, filesWrote[bestRunIndex]+".fasta.coords");
    float[][] Solution2 = new float[numNodes][2];
    float[][] rotateMemory = new float[numNodes][2];
    ana.println("Start of rotations:");
    for(int run = 0; run < numTimes; run++){
      if (run==bestRunIndex) continue; //Very necessary clause!!!
      getCoordsFile(Solution2, filesWrote[run]+".fasta.coords");
      //So, modify filesWrote[run] so that it's better.
      ana.printf("%10d",run);
      rotateToBest(Solution1, Solution2, rotateMemory, filesWrote[run]+".fasta", ana);
    }
    ana.close();
  } 
  catch (Throwable e){
    e.printStackTrace();
    batchError("Error in analysis:"+e.toString());
  }
  thenExit();
}
private void thenExit(){
  System.exit(0);
}
private void batchError(String msg){
  System.out.println();
  System.out.println("Try -Dhelp to read the help file!");
  System.out.println("//////////////////////////////////PROBLEM");
  System.out.println("|>   "+msg);
  System.out.println("//////////////////////////////////");
  System.out.println();
  System.out.println("I Encountered the problem when...");
  throw new RuntimeException("Batch error.");
}
private void doHelp(){
  System.out.println("////////////////////////////////");
  System.out.println("|>      Dgraph batch-file help:");
  System.out.println("////////////////////////////////");
  String formatString = "%-10s     %s";
  System.out.printf(formatString,"-Druns","How many times to run the simulation");
  System.out.println();
  System.out.printf(formatString,"-Diters","How many iterations each run involves");
  System.out.println();
  System.out.printf(formatString,"-Dout","A directory in which to save the results");
  System.out.println();
  System.out.printf(formatString,"-Dprefix","A filename prefix for all result files");
  System.out.println();
  System.out.printf(formatString,"-Dparams","6 numeric simulation parameters seperated by commas.\r\n\t *In the order they appear when you press 'j' in the GUI.");
  System.out.println();
  System.out.printf(formatString,"-Din","The fasta file to read");
  System.out.println();
  System.out.printf(formatString,"-DinCustom","An NxN matrix (in a text file) of distances to use");
  System.out.println();
  System.out.printf(formatString,"-Dpdf","Flag, means to open the fasta, render a pdf in the same directory, and exit.");
  System.out.println();
  System.out.printf(formatString,"-Dzoom","Zooming factor. Lower numbers means more zoomed in. Must be >=1.");
  System.out.println();
}










