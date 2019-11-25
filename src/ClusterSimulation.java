import java.awt.Color;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeSet;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextArea;

import processing.core.PApplet;
import processing.core.PFont; 

public class ClusterSimulation extends PApplet {
String removed = "";

boolean PRESETUP = false;
int WHAT_DO_YOU_WANT = -1;
boolean singleFramePdf = false;

public void setup() {
	  doCommandLine();
	  if (singleFramePdf){
	    PRESETUP = true;
	    WHICH_SCREEN = 2;
	  }
	  if (!PRESETUP){
	    PRESETUP = true;
	    //Show a dialog with options of what to do:
	    String [] optionsStr = { 
	      "Run DGraph", 
	      "Utility: Sequence list -> Fasta Format", 
	      "Utility: Fasta Format -> Sequence List", 
	      "Utility: Remove Duplicates from Sequence List", 
	      "Utility: Fasta Format -> List of 'Fasta Headers'",
	      "Utility: Clustalw output -> Alignment score distance matrix",
	      "Utility: Sequence List -> Neo-PD Score Matrix (Long Proteins)",
	      "Utility: Sequence List -> PD Score Matrix (Short Peptide)",
	      "Utility: DNA Fasta -> DNA Sequence Similarity matrix",
	      "Utility: Find all pairs of peptide in a list with PD Score under a threshold",
	      };
	    JRadioButton [] options = new JRadioButton[optionsStr.length];
	    ButtonGroup group = new ButtonGroup();
	    JPanel pane = new JPanel();
	    pane.setLayout(new GridLayout(0,1));
	    for(int k = 0; k < optionsStr.length;k++){
	      options[k] = new JRadioButton (optionsStr[k]);
	      group.add(options[k]);
	      pane.add(options[k]);
	    }
	    options[0].setSelected(true);

	    Frame holder = new Frame("DGraph - Initializing.");
	    holder.setSize(200,0);
	    holder.setVisible(true);
	    holder.requestFocus();
	    try {
	      Thread.sleep(100);
	    } 
	    catch (Throwable e){
	    };
	    holder.requestFocus();

	    boolean shouldRun = JOptionPane.showOptionDialog(
	    holder,
	    pane,
	    "Welcome! What would you like to run?",
	    JOptionPane.OK_CANCEL_OPTION,
	    JOptionPane.PLAIN_MESSAGE,null,null,null) == JOptionPane.OK_OPTION;
	    if (!shouldRun) {
	    	System.exit(0);
	    }

	    holder.setVisible(false);

	    int choice = 0;
	    for(int k = 0; k < options.length; k++){
	      if (options[k].isSelected()){
	        choice = k;
	        break;
	      }
	    }
	    WHAT_DO_YOU_WANT = choice-1; //0=-1
	  }
	  
	  //Test here if we have write priviledges
	  String test = sketchPath+File.separator+".test";
	  try {
	    FileOutputStream fos = new FileOutputStream(test);
	    fos.write(2);
	    fos.close();
	  } 
	  catch (Throwable e){
	    JOptionPane.showMessageDialog(frame,"Please copy this program to a directory that has file-write access.\nNo installation is necessary except for this step.\n"+/*sketchPath+*/" does not have write access.","Write-Protected Directory",JOptionPane.ERROR_MESSAGE);
	    System.exit(1);
	  }
	 
	  if (WHAT_DO_YOU_WANT==-1){
	    size(800,800,P2D);
	  } 
	  else {
	    TestFasta(WHAT_DO_YOU_WANT);
	    System.exit(0);
	  }
    txt = createFont("Vera.ttf",90,true);
    //textPaneClone = createGraphics(2048,2048,P2D);
    //textPaneClone = createGraphics(width,height,P2D);
    //beginDraw - endDraw need to be tight.
    //textPaneClone.beginDraw();
    
    frameRate(60);
}

int WHICH_SCREEN = 0;
public void draw(){
    if (!frame.isResizable()){
      frame.setResizable(true);
    }
	if (WHICH_SCREEN == 1){
		drawSimulation(true);
		if (keyPressed && key=='p'){
			destroyShaded();
			new Thread() {
				public void run() {
					setupShaded();
				}
			}.start();
			WHICH_SCREEN++;
		}
	} else if (WHICH_SCREEN==2){
		drawShaded();
	} else {
		background(0,0,100);
		destroySimulation();
		new Thread() {
			public void run() {
				setupSimulation();
			}
		}.start();
		WHICH_SCREEN = 1;
	}
}

float[] shadings;
int[] nodeConnections;
long gotHere;
String pdfUrl = null;
int shadedScreen = 0;
//Output directory for drawShaded
String outputPdfDir = null;
public void destroyShaded() {
	shadedScreen = 0;
}
public void setupShaded() {
	gotHere = System.nanoTime();
	shadings = new float[nodes.length]; //defaults to 0
	nodeConnections = new int[nodes.length];

	String toRead = null;
	String[] commands = null;
	if (singleFramePdf){
		commands = new String[0];
	}
	if (true) throw new RuntimeException();
	if (commands==null){
		FileDialog fd = new FileDialog(frame, "Open the coloring file", FileDialog.LOAD);
		//while ( fd.getDirectory() == null && fd.getFile() == null) 
		//{
		fd.show();
		//}
		try {
			outputPdfDir = fd.getDirectory().toString();
			File got = new File( fd.getDirectory() + File.separator + fd.getFile());
			toRead = got.getAbsoluteFile().toString();
		} 
		catch (Throwable e){
			JOptionPane.showMessageDialog(frame,"Incorrect file format");
			e.printStackTrace();
		};
		commands = loadStrings(toRead);
	}
	for(int k = 0; k < commands.length; k++){
		String[] line = commands[k].trim().split("\\s+");
		if (line.length<2) {
			continue;
		}
		for(int p = 0; p < nodes.length; p++){
			if (nodes[p].FastaLabel.equalsIgnoreCase(line[0])){
				shadings[p] = PApplet.parseFloat(line[1]);
				if (line.length>=3){
					for(int q = 0; q < nodes.length; q++){
						if (nodes[q].FastaLabel.equalsIgnoreCase(line[2])){
							nodeConnections[p] = q;
						}
					}
				} else {
					nodeConnections[p] = -1;
				}
				break;
			}
		}
	}
	//This fires off drawShaded to start doing work.
	shadedScreen = 1;
}
public void drawShaded(){
	if (shadedScreen == 0) {
		return;
	} 
	boolean renderThisFrame = shadedScreen == 1;
	shadedScreen = 2;

	if (renderThisFrame){
		beginRecord(PDF, pdfUrl!=null?pdfUrl:(outputPdfDir+"/"+frameCount+".pdf"));
	}
	for(int k = 0; k < nodes.length; k++){
		if (!nodes[k].exists) continue;
		nodes[k].shadeColor = lerpColor(color(0,0,0),color(255,0,0),constrain(shadings[k],0,1)); //erased after 1 use, by the way.
		nodes[k].nodeConnection = nodeConnections[k];
	}
	//Superdraw:
	boolean eatKeyPress = keyPressed;
	keyPressed = false;
	drawSimulation(false);
	keyPressed = eatKeyPress;
	if (renderThisFrame){
		endRecord();
		if (singleFramePdf){
			System.exit(0);
		}
	}
	textFont(txt);
	textSize(txtFontSizeMedium * txtFontSizeScaling);
	fill(0);
	textAlign(LEFT,TOP);
	text("Press any key to return.",0,0);

	if (keyPressed && (System.nanoTime()-gotHere)>.5e9f){
		WHICH_SCREEN = 1;
		shadings = null;
	}
}
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
public void doCommandLine(){
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
  setupSimulation();

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



public static final int SHOW_LABELS_SEQUENCE = 0, 
SHOW_LABELS_NAME = SHOW_LABELS_SEQUENCE + 1, 
SHOW_LABELS_NUMBER = SHOW_LABELS_NAME + 1, 
SHOW_LABELS_NONE = SHOW_LABELS_NUMBER + 1;

private String toRead = null;
private String scoreName = "";
public String pickUserFile(String title) {
  FileDialog fd = new FileDialog(frame, title, FileDialog.LOAD);
  fd.show();
  if ( fd.getDirectory() != null && fd.getFile() != null) 
  {
    try {
      File got = new File( fd.getDirectory() + File.separator + fd.getFile());
      if (!got.exists()) {
        JOptionPane.showMessageDialog(frame, "You entered a nonexistant file:\n"+got);
        return pickUserFile(title);
      }
      return got.getAbsoluteFile().toString();
    } 
    catch (Throwable e) {
      e.printStackTrace();
    };
    //println("You chose to open this file: " + chooser.getSelectedFile().getName());
  } 
  return "data.txt";
}
float[][] distanceMatrix = null; //If null, calc PD values.    
public boolean iS = false; //iS means inScript. So, we have to ignore all GUI updaters in ReInit.
public void destroySimulation() {
	  simulationTicks = 0;
}
public void setupSimulation() {
  //if (!iS) g.endDraw();

  /*
  if (toRead!=null && !(new File(toRead).exists())){
   toRead = null;
   }
   */
  if (toRead==null) {
    toRead = pickUserFile("Open a Fasta File");
    boolean hasCustomMatrix = JOptionPane.showConfirmDialog(frame, "Do you want to use a precalculated distance matrix?\n(Choose no and I will calculate the pairwise PD-distances for you!)", "Custom Distances", JOptionPane.YES_NO_OPTION)==JOptionPane.YES_OPTION;

    if (hasCustomMatrix) {
      try {
        String customLoc = pickUserFile("Input a Custom Matrix");
        GetCustomDistances(customLoc);
      } 
      catch (ArrayIndexOutOfBoundsException e) {
        JOptionPane.showMessageDialog(frame, "Matrix was not square.");
        exit();
      }
    }
  }

  //if (!iS) g.beginDraw();

  //System.out.println(toRead);
  String[] lines = loadStrings(toRead);
  String[][] fastaInfo = null;
  try {
    fastaInfo = unFastaSuper(lines);
  } 
  catch (Throwable e) {
    JOptionPane.showMessageDialog(frame, "Invalid Fasta File.", "Error!", JOptionPane.ERROR_MESSAGE);
    System.exit(1);
  }
  Node[] nodes = new Node[fastaInfo[0].length];
  removed = "";
  
  //if (!iS) noLoop();
  
  long lastDraw = System.nanoTime();
  if (distanceMatrix==null) {
    scoreName = "PD";
  } 
  else {
    scoreName = JOptionPane.showInputDialog(frame, "Enter the name of your custom score function");
  }
  for (int k = 0; k < nodes.length; k++) {
    nodes[k] = new Node();
  }
  this.nodes = nodes;
  for (int k = 0; k < nodes.length; k++) {
    nodes[k].FastaLabel = fastaInfo[0][k];
    nodes[k].FastaSequence = fastaInfo[1][k];
    final String hashCoordText = "#COORDS:";
    int hashCoords = nodes[k].FastaLabel.indexOf(hashCoordText);
    if (hashCoords==-1) { //Random position
      scatterNode(nodes[k]);
    } 
    else { //Nonrandom position, strip the coords for the actual header.
      String[] pos = nodes[k].FastaLabel.substring(hashCoords+hashCoordText.length()).split(","); //Split on comma for x,y
      nodes[k].pos[0] = new Float(pos[0]);
      nodes[k].pos[1] = new Float(pos[1]); 
      float PIX_TO_PD_FACTOR = (ZOOM/640.f);
      nodes[k].pos[0]/=PIX_TO_PD_FACTOR;
      nodes[k].pos[1]/=PIX_TO_PD_FACTOR;
      nodes[k].FastaLabel = nodes[k].FastaLabel.substring(0, hashCoords).trim(); //The "Real" header, remove trailing leading whitespace.
      //System.out.printf("Initial position: %s (%.3f,%.3f)",nodes[k].FastaLabel,nodes[k].pos[0],nodes[k].pos[1]);
      //System.out.println();
    }
    float[] dist = new float[nodes.length];
    float dontShowCutoff = REMOVAL_ISLAND_FINDER;
    boolean hasLessThanCutoff = false;
    for (int p = 0; p < nodes.length; p++) {
      if (p==k) continue;
      if (distanceMatrix==null) {
        dist[p]=PDScore(fastaInfo[1][p], fastaInfo[1][k]);
      } 
      else {
        if (distanceMatrix.length < nodes.length) {  
          JOptionPane.showMessageDialog(frame, "Distance Matrix did not match Fasta File in length"+"("+distanceMatrix.length+" "+nodes.length+")");
          exit();
        }
        dist[p]=distanceMatrix[k][p];
      }
      if (dist[p] < dontShowCutoff /**&& dist[p]!=0**/)
        hasLessThanCutoff = true;
    }
    nodes[k].distances = dist;
    nodes[k].myId = k;
    nodes[k].exists = hasLessThanCutoff;
    if (!hasLessThanCutoff) {
      removed+=" ,"+(k+1);
    }
    //System.out.println("Processed: "+lines[k]);
    if (!iS) if ((System.nanoTime()-lastDraw)>.5e9f) {
      lastDraw = System.nanoTime();
      background(0);
      fill(255);
      textFont(txt);
      textAlign(LEFT, TOP);
      textSize(txtFontSizeMedium * txtFontSizeScaling);
      text("Calculating PDscores:"+nf((k/(float)nodes.length)*100, 0, 2)+"%", 5, 5);
      g.endDraw();
      //repaint();
      g.beginDraw();
    }
  }
  //Phew.
  
  if (!iS) {
    //loop();
    //g.endDraw();
    /*
    frame.setSize(new Dimension(frame.getWidth()-10, frame.getHeight()-10));
    try {
      Thread.sleep(100);
    }
    catch (Throwable e) {
    };
    frame.setSize(new Dimension(frame.getWidth()+10, frame.getHeight()+10));
    */
    //g.beginDraw();
  }
  
  simulationTicks = 1;
}

public void scatterNode(Node got) {
  float theta = random(TWO_PI);
  float mag = random(640/200f, 640/2f);
  got.pos[0] = mag*cos(theta);
  got.pos[1] = mag*sin(theta);
}
public void scatterNodes() {
  //Similar to reinit, except that only the positions of the already existant nodes are moved.
  for (int k = 0; k < nodes.length; k++) {
    Node got = nodes[k];
    if (!got.exists) continue;
    scatterNode(got);
  }
}
public void strokeColors(float val) {
  float var = constrain(val/SHOW_LINES, 0, 1);
  float amtRed = pow(constrain(-2*var+1+.4f, 0, 1), 2);
  float amtYello = constrain(1-3*abs(var-.3f), 0, 1);
  float amtBlue = sqrt(constrain(-1.5f*(1-var)+1, 0, 1));
  int c = lerpColor(lerpColor(color(0, 0, 255), color(0, 170, 40), (1-amtRed)), color(255, 255, 152), amtBlue);
  stroke(lerpColor(c, color(255, 255, 255), var*.5f));
}
private Node[] nodes;
private float TOTAL_VARIANCE;
private float RANDOM_VARIANCE = -1;
private final float ZOOM = 50; //ZOOM effects the scaling of the position integration; it used to be determined by screen coords ...
private float VIEW_SCALE = -1; //Updated from VIEW_SCALE_USER
private float VIEW_SCALE_USER = 25;
private int showLabels = SHOW_LABELS_NONE; 
private boolean ShowMinLines = false;
private int NodeDotSize = 3; //rectangle of this radius
private String underMouse;
int simulationTicks;

private boolean UPDATES_ON_DRAW = true; //Normally true, but set to false on the first frame in the GUI.
private boolean wantsFrameBreak = false;
private boolean simStarted = false;
public void drawSimulation(boolean controlInfo) {
  if (simulationTicks == 0) {
	  return;
  }
  UPDATES_ON_DRAW = simStarted; //GUI uses this so the user can just look at it for a second.
  underMouse = "";
  background(255);
  textFont(txt);
  textSize(txtFontSizeMedium * txtFontSizeScaling);
  noStroke();
  fill(0);
  if (controlInfo) {
    textAlign(LEFT, TOP);
    text("\"Islands\" Removed by "+REMOVAL_ISLAND_FINDER+" PD threshold: "+removed.substring(min(removed.length(), 2)), 2, controlInfo?20:5, width, 80);
    text("Noncomparable Seperation Factor:"+CURRENT_CUTOFF, 2, 4);
    text("Displaying:"+SHOW_LINES, 2, 100);
    textAlign(RIGHT, TOP);
    if (!UPDATES_ON_DRAW) {
      fill(255, 0, 0);
      text("Simulation stopped. Press 'y' to start again.", width, 0);
    }
    fill(0);
    int fromBottomRight = 150;
    text("'h' to switch on/off lines", width-5, height-fromBottomRight-33);
    text("'c' to turn off labels", width-5, height-fromBottomRight-20);
    text("'s' to show sequence numbering", width-5, height-fromBottomRight-7);
    text("'a' to show sequences", width-5, height-fromBottomRight+5);
    text("'b' to show sequence headers", width-5, height-fromBottomRight+18);
    text("'9' zooms out", width-5, height-fromBottomRight+35);
    text("'0' zooms in", width-5, height-fromBottomRight+45);
    text("'m' saves the current state", width-5, height-fromBottomRight+55);
    text("'p' proceeds to the analysis step", width-5, height-fromBottomRight+72);
    text("'j' sets the physics options", width-5, height-fromBottomRight+85);
    text("'d' disable nodes", width-5, height-fromBottomRight+100);
    text("'w' increases dot-size", width-5, height-fromBottomRight+115);
    text("'q' decreases dot-size", width-5, height-fromBottomRight+125);
    text("'r' or 'e' re-run the simulation. 'e' only scatters the points, 'r' does more.", width-5, height-fromBottomRight+135);
  }
  if (true) {
    if (transparentCountdown>0) {
      nodes[0].setStroke(8, 0);
      text("Overstretched by 8", 0, 150);
      line(0, 150, 50, 150);
      nodes[0].setStroke(0, 0);
      text("Just about perfect", 0, 175);
      line(0, 175, 50, 175);
      nodes[0].setStroke(-8, 0);
      line(0, 200, 50, 200);
      text("Understretched by 8", 0, 200);
    } 
    else {
      textAlign(LEFT, CENTER);
      text(scoreName+" SCORE - DISTANCE", 20, 130);
      beginShape(QUADS);
      float minY = 150;
      float maxY = 400;
      for (int PD = 0; PD <= SHOW_LINES; PD++) {
        float mag = PD/SHOW_LINES;
        float y = lerp(minY, maxY, mag);
        strokeColors(PD);
        fill(g.strokeColor);
        if (PD>0) {
          vertex(20, y);
          vertex(25, y);
        }
        vertex(25, y);
        vertex(20, y);
      }
      endShape();
      fill(0);
      for (int PD = 0; PD <= SHOW_LINES; PD++) {
        float mag = PD/SHOW_LINES;
        float y = lerp(minY, maxY, mag);
        text(PD, 30, y);
      }
    }
  }
  //Draw simulation:
  VIEW_SCALE = 40f/VIEW_SCALE_USER;
  //keyPressed2();
  //pushMatrix();
  translate(width/2, height/2);
  //The actual simulation:
  runActualSimulation();
  //End actual simulation.
  flushLines();
  if (RANDOM_VARIANCE<0) {
    RANDOM_VARIANCE = TOTAL_VARIANCE;
  }

  labelPlanner.clear();
  for (int k = 0; k < nodes.length; k++) {
	  nodes[k].drawDot();
  }
  for (int k = 0; k < nodes.length; k++) {
	  nodes[k].draw();
  }
  wantsFrameBreak = false;
  //popMatrix();
  transparentCountdown--;
  fill(0);
  textSize(txtFontSizeMedium * txtFontSizeScaling);
  textAlign(LEFT, BOTTOM);
  translate(-width/2, -height/2);
  if (controlInfo) {
    text("Simulation ticks:"+simulationTicks+"\nSolution \'goodness\':"+TOTAL_VARIANCE+"\nTimes smaller than random:"+RANDOM_VARIANCE/TOTAL_VARIANCE+"\nNodes under mouse:"+underMouse.substring(min(underMouse.length(), 2)), 0, height-200, min(300, width-100), 200);
    stroke(0);
    if (false && simulationTicks==200) {
      System.out.printf("%-7.3f %-7.3f %-7.3f %-14.3f", timeStep, fluidFriction, mass, TOTAL_VARIANCE);
      System.out.println();
    }
  }
  /*
  if (frameCount % 2 == 0) {
	  image(textPaneClone,0,0,width,height,0,0,width,height);
  }
  */
}
private long acted = System.nanoTime();
/**
 * Wrapper around _keyPressed that throttles actions and interrupts ongoing draw()
 */
public void keyPressed() {    
  wantsFrameBreak = true;
  //key = e.getKeyChar();
  keyPressed = true;
  if (System.nanoTime() - acted > .033e9f) {
	  _keyPressed();
  }
}
public void _keyPressed() {
  if (keyPressed && Character.toLowerCase(key)=='y') {
    acted = System.nanoTime();
    simStarted = !simStarted;
    return;
  }
  if (keyPressed && Character.toLowerCase(key)=='d') {
    acted = System.nanoTime();
    boolean canDo = JOptionPane.showConfirmDialog(frame, "Disable all nodes under your mouse? This cannot be undone.")==JOptionPane.OK_OPTION;
    if (!canDo) {
    	return;
    }
    String[] nums = underMouse.split("[, ]+");
    for (int k = 0; k < nums.length; k++) {
      try {
        int o = Integer.parseInt(nums[k].trim())-1;
        nodes[o].exists = false;
      } 
      catch (Throwable e) {
      }
    }
    return;
  }
  if (keyPressed && Character.toLowerCase(key)=='r') {
    acted = System.nanoTime();
    setupSimulation();
    return;
  }
  if (keyPressed && Character.toLowerCase(key)=='e') {
    acted = System.nanoTime();
    scatterNodes();
    return;
  }
  if (keyPressed && Character.toLowerCase(key)=='m') {
    acted = System.nanoTime();
    saveNodeStates(null); //Prompt dialog
    return;
  }
  if ((keyPressed && Character.toLowerCase(key)=='t')) {
    acted = System.nanoTime();
    transparentCountdown = 10;
    return;
  }
  /*
  if (frame!=null && keyPressed && Character.toLowerCase(key)=='c'){
   acted = System.nanoTime();
   save("Screenshots/"+(int)SHOW_LINES+"-"+frameCount+".png");
   SHOW_LINES--;
   try {
   Thread.sleep(200);
   } 
   catch (Throwable e){
   };
   return;
   }
   */

  if ((keyPressed && Character.toLowerCase(key)=='a')) {
    acted = System.nanoTime();
    showLabels = SHOW_LABELS_SEQUENCE;
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='b')) {
    acted = System.nanoTime();
    showLabels = SHOW_LABELS_NAME;
    return;
  }


  if ((keyPressed && Character.toLowerCase(key)=='s')) {
    acted = System.nanoTime();
    showLabels = SHOW_LABELS_NUMBER;
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='c')) {
    acted = System.nanoTime();
    showLabels = SHOW_LABELS_NONE;
    return;
  }


  if ((keyPressed && Character.toLowerCase(key)=='h')) {
    acted = System.nanoTime();
    ShowMinLines = !ShowMinLines;
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='w')) {
    acted = System.nanoTime();
    NodeDotSize=min(NodeDotSize+1, 10);
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='q')) {
    acted = System.nanoTime();
    NodeDotSize=max(NodeDotSize-1, 1);
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='9')) {
    acted = System.nanoTime();
    VIEW_SCALE_USER=max(VIEW_SCALE_USER-1, 1);
    return;
  }


  if ((keyPressed && Character.toLowerCase(key)=='0')) {
    acted = System.nanoTime();
    VIEW_SCALE_USER=VIEW_SCALE_USER+1;
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='7')) {
    acted = System.nanoTime();
    txtFontSizeScaling=max(txtFontSizeScaling - 0.05f,0.25f);
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='8')) {
    acted = System.nanoTime();
    txtFontSizeScaling=min(txtFontSizeScaling + 0.05f,4f);
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='j')) {
    acted = System.nanoTime();
    newPhysicsVariables();
    return;
  }
}
public void runActualSimulation() {
  TOTAL_VARIANCE = 0; //Calculate
  if (UPDATES_ON_DRAW) {
    simulationTicks++;
  }
  if ((lineRendering==null || lineRendering.length < nodes.length*nodes.length) && nodes.length < 3000) {
    lineRendering = new LineSort[nodes.length*nodes.length];
    for (int k = 0; k < lineRendering.length; k++) {
      lineRendering[k] = new LineSort();
    }
  }
  for (int k = 0; k < nodes.length; k++) {
    nodes[k].update();
  }
  float[] averagePos = new float[2];
  int realNodes = 0;
  for (int k = 0; k < nodes.length; k++) {
    if (!nodes[k].exists) continue;
    realNodes++;
    averagePos[0] += nodes[k].pos[0];
    averagePos[1] += nodes[k].pos[1];
  }
  averagePos[0]/=realNodes;
  averagePos[1]/=realNodes;
  for (int k = 0; k < nodes.length; k++) {
    nodes[k].pos[0] -= averagePos[0];
    nodes[k].pos[1] -= averagePos[1];
  }
}
public class LineSort implements Comparable {
  public float myz;
  float[] myLine = new float[4];
  int strokeColor;
  boolean dead = true;
  public int compareTo(Object oth) {
    return (int)Math.signum(myz-((LineSort)oth).myz);
  }
}
int lineRenderingUsed = 0;
private LineSort[] lineRendering;
public void flushLines() {
  strokeWeight(1);
  if (lineRendering==null) return;
  Arrays.sort(lineRendering, 0, lineRenderingUsed);
  for (int k = 0; k < lineRenderingUsed; k++) {
    LineSort ls = lineRendering[k];
    if (ls.dead) break;
    ls.dead = true;
    float[] dat = ls.myLine;
    stroke(ls.strokeColor);
    line(dat[0], dat[1], dat[2], dat[3]);
  }

  lineRenderingUsed = 0;
}
int transparentCountdown = 0;
private int selected = -1;
PFont txt;
//In pt
final float txtFontSizeSmall = 12;
final float txtFontSizeMedium = 20;
float txtFontSizeScaling = 1.0f;

/** DISTANCE VARIABLES **/

private float CURRENT_CUTOFF = 999;
private float SHOW_LINES = 14;
private float REMOVAL_ISLAND_FINDER = 12;

/** PHYSICS VARIABLES **/

float timeStep = .001f;
float fluidFriction = .008f;
float mass = .01f;

abstract class VarSetter extends JPanel {
  public VarSetter(int numSettings) {
    this.numSettings = numSettings;
    areas = new JTextArea[numSettings];
    labels = new JLabel[numSettings];
  }
  int numSettings;
  JTextArea[] areas;
  JLabel[] labels;
  public abstract void dedicate();
}
class PhysicsSetter extends VarSetter {
  String[] names = new String[] {
    "Delta T", "Fluid Friction", "Mass"
  };
  String[] initValues = new String[] {
    timeStep+"", fluidFriction+"", mass+""
  };
  public PhysicsSetter() {
    super(3);
    setLayout(new GridLayout(numSettings, 1));
    for (int k = 0; k < areas.length; k++) {
      areas[k] = new JTextArea();
      areas[k].setText(initValues[k]);
      labels[k] = new JLabel();
      labels[k].setLabelFor(areas[k]);
      labels[k].setText(names[k]);
      labels[k].setVerticalTextPosition(JLabel.BOTTOM);
      labels[k].setHorizontalTextPosition(JLabel.LEFT);
      areas[k].setBorder(BorderFactory.createLineBorder(Color.black));
      add(labels[k]);
      add(areas[k]);
    }
  }
  public void dedicate() {
    timeStep = PApplet.parseFloat(areas[0].getText());
    fluidFriction = constrain(PApplet.parseFloat(areas[1].getText()), 0, 1);
    mass = PApplet.parseFloat(areas[2].getText());
  }
}
class DistanceSetter extends VarSetter {
  String[] names = new String[] {
    "Comparison Cutoff", "Display Lines < than", "Island Removal Cutoff"
  };
  String[] initValues = new String[] {
    CURRENT_CUTOFF+"", SHOW_LINES+"", REMOVAL_ISLAND_FINDER+""
  };
  public DistanceSetter() {
    super(3);
    setLayout(new GridLayout(numSettings, 1));
    for (int k = 0; k < areas.length; k++) {
      areas[k] = new JTextArea();
      areas[k].setText(initValues[k]);
      labels[k] = new JLabel();
      labels[k].setLabelFor(areas[k]);
      labels[k].setText(names[k]);
      labels[k].setVerticalTextPosition(JLabel.BOTTOM);
      labels[k].setHorizontalTextPosition(JLabel.LEFT);
      areas[k].setBorder(BorderFactory.createLineBorder(Color.black));
      add(labels[k]);
      add(areas[k]);
    }
  }
  public void dedicate() {
    CURRENT_CUTOFF = PApplet.parseFloat(areas[0].getText());
    SHOW_LINES = PApplet.parseFloat(areas[1].getText());
    REMOVAL_ISLAND_FINDER = PApplet.parseFloat(areas[2].getText());
  }
}
public void newPhysicsVariables() {
  VarSetter options1 = new PhysicsSetter();
  VarSetter options2 = new DistanceSetter();
  JPanel options = new JPanel();
  options.setLayout(new GridLayout(1, 2));
  options.add(options2);
  options.add(options1);
  boolean changedSettings = JOptionPane.showConfirmDialog(frame, options, "Settings", JOptionPane.OK_CANCEL_OPTION)==JOptionPane.OK_OPTION;
  if (changedSettings) {
    options1.dedicate();
    options2.dedicate();
  }
}

/** END PHYSICS **/

class LabelPlanner {
	private boolean[][] markedPoints;
	public void clear() {
		if (markedPoints == null || markedPoints.length != height * width) {
			markedPoints = new boolean[height][width];
		}
		for(boolean[] x : markedPoints) {
			Arrays.fill(x, false);
		}
	}
	public boolean canDrawAt(float ox, float oy) {
		if (ox < 0 || ox >= markedPoints[0].length || oy < 0 || oy >= markedPoints.length) return false;
		boolean toRet = !markedPoints[(int)oy][(int)ox];
		return toRet;
	}
	public void mark(float ox, float oy, float r) {
		//Include 1 pixel of border:
		r+=1;
		
		int maxX = (int)(ox + r + 1);
		int maxY = (int)(oy + r + 1);
		float rs = r * r;
		for(int y = (int)(oy - r); y < maxY; y++) {
			float dys = (y - oy)*(y - oy);
			if (y < 0 || y >= markedPoints.length) continue;
			for(int x = (int)(ox - r); x < maxX; x++) {
				if (x < 0 || x >= markedPoints[0].length) continue;
				float dxs = (x - ox)*(x - ox);
				if (dxs + dys < rs) {
					markedPoints[y][x] = true;
				}
			}
		}
	}
}
LabelPlanner labelPlanner = new LabelPlanner();

class Node {
  float acc [] = new float[2];
  float vel [] = new float[2];
  float pos [] = new float[2];
  float[] distances;
  boolean exists;
  int shadeColor = -1;
  int nodeConnection = -1;
  String FastaLabel;
  String FastaSequence;
  int myId;
  public void draw() {
    if (!exists) {
      return;
    }
    float x = pos[0]*VIEW_SCALE;
    float y = pos[1]*VIEW_SCALE;
    switch (showLabels) {
    case SHOW_LABELS_SEQUENCE:
      placeText(FastaSequence, x, y);
      break;
    case SHOW_LABELS_NAME:
      placeText(FastaLabel, x, y);
      break;
    case SHOW_LABELS_NUMBER:
      placeText(myId+1+"", x, y);
    }
  }
  public void drawDot() {
    if (!exists) {
      return;
    }
    float x = pos[0]*VIEW_SCALE;
    float y = pos[1]*VIEW_SCALE;
    strokeWeight(1);
    if (shadeColor!=-1) {
      fill(shadeColor);
      stroke(shadeColor);
      shadeColor = -1;
    } 
    else {
      stroke(0);
      fill(255);
    }
    if (nodeConnection!=-1) {
      drawLineTo0(nodeConnection, brightness(shadeColor)/60.f);
      nodeConnection = -1;
    }
    rect(x-NodeDotSize, y-NodeDotSize, NodeDotSize*2, NodeDotSize*2);
    float cx = screenX(0, 0);
    float cy = screenY(0, 0);
    labelPlanner.mark(cx + x, cy + y, NodeDotSize);
  }
  public void placeText(String labelString, float offx, float offy) {
    float cx=0, cy=0;
    
    fill(0);
    textFont(txt);
    textSize(txtFontSizeSmall * txtFontSizeScaling);
    textAlign(CENTER, CENTER);
    translate(offx, offy);
      
    
    cx = screenX(0, 0);
    cy = screenY(0, 0);
    if (pow(mouseX - cx, 2)+pow(mouseY-cy, 2) < 3000) {
      translate(-offx, -offy);
      underMouse+=", "+(myId+1);
      return;
    }

    //Default to an angle that's away from the center of the picture (which is at 0,0)
    float angle0 = atan2(offy, offx);
    //Try multiple angles for the label ray until we find one that doesn't conflict with other labels
    float angleDelta = TWO_PI / 10;
    //Label will be placed distance r away, at some angle, in a circle of radius circleR
    float circleR = 20;
    float r = max(25, NodeDotSize + 5);
    for (float angle = angle0; angle < angle0+TWO_PI; angle += angleDelta) {
      boolean canDraw = true;
      
      //Angles close to horizontal are ugly, and lead to unreadable labels when showing sequence
      //as these overflow the "circleR" space at the label. So avoid that:
      float horizThreshold = TWO_PI/20;
      if (abs(angle + horizThreshold) % PI < horizThreshold * 2) {
        canDraw = false;
      }
      
      //Check drawability slightly to the left and right of the proposed label position
      if (!labelPlanner.canDrawAt(cx + r*cos(angle) - circleR, cy + r*sin(angle))) {
          canDraw = false;
      }
      if (!labelPlanner.canDrawAt(cx + r*cos(angle) + circleR, cy + r*sin(angle))) {
          canDraw = false;
      }
      //Check drawability on the segment connecting label to point. Can't check at the center of the ray since we want to allow
      //labels for points radiating from the exact same spot.
      for (float rad = max(10, NodeDotSize + 5); rad < r + circleR && canDraw; rad+= 5) {
        float ox = rad*cos(angle);
        float oy = rad*sin(angle);
        if (!labelPlanner.canDrawAt(cx + ox, cy + oy)) {
            canDraw = false;
        }
      }
      if (canDraw) {
        float ox = r*cos(angle);
        float oy = r*sin(angle);
        
        stroke(0, 0, 0, 100);
        line(0, 0, ox*((r-circleR/2)/r), oy*((r-circleR/2)/r));
        fill(1, 1, 1, 40);
        if (showLabels == SHOW_LABELS_NUMBER) {
          ellipse(ox, oy, circleR, circleR);
        }
        fill(0);
        translate(-offx, -offy);
        text(labelString, (int)(offx+ox), (int)(offy+oy));
        labelPlanner.mark(cx + ox, cy + oy, circleR);
        return;
      }
    }
    
    translate(-offx, -offy);
  }
  public void setStroke(float goodness, float distance) {
    strokeWeight(2);
    colorMode(HSB, 255);
    stroke(200+goodness/SHOW_LINES*50, 255, 255);
    colorMode(RGB, 255);
  }
  public boolean drawsLines() {
    return (transparentCountdown>0 || ShowMinLines) && !wantsFrameBreak;
  }
  public void drawLineTo(int oth, float goodness) {
    //called by update, but if we're in script, don't do it!
    if (iS) return;
    //Done.
    boolean trans = transparentCountdown>0;
    if (!trans) {
      if (ShowMinLines && distances[oth]<SHOW_LINES) {
        strokeWeight(2);
        strokeColors(distances[oth]);
      } 
      else {
        return;
      }
    } 
    else {
      setStroke(goodness, distances[oth]);
    }

    drawLineTo0(oth, goodness);
  }
  public void drawLineTo0(int oth, float goodness) {
    float x = pos[0]*VIEW_SCALE;
    float y = pos[1]*VIEW_SCALE;
    Node other = nodes[oth];
    float ax = other.pos[0]*VIEW_SCALE;
    float ay = other.pos[1]*VIEW_SCALE;
    float z = -distances[oth]/(float)SHOW_LINES*8;


    if (lineRendering!=null && lineRenderingUsed<lineRendering.length) {
      LineSort ls = lineRendering[lineRenderingUsed];
      ls.strokeColor = g.strokeColor;
      ls.myLine[0]=x;
      ls.myLine[1]=y;
      ls.myLine[2]=ax;
      ls.myLine[3]=ay;
      ls.myz = z;
      ls.dead = false;
      lineRenderingUsed++;
    }

    /*
    line(x,y,ax,ay);
     flush();
     */
    colorMode(RGB, 255);
  }
  public void update() {
    if (!exists) {
      return;
    }

    if (!mousePressed) {
      selected = -1;
    }
    if (mousePressed) {
      if (selected==myId) {
        float[] disp = componentsTo(new float[] {
          (mouseX-width/2)/VIEW_SCALE, (mouseY-height/2)/VIEW_SCALE
        }
        );
        addVec(pos, disp);
        selected = -1;
      }
      float x = pos[0]*VIEW_SCALE+width/2;
      float y = pos[1]*VIEW_SCALE+height/2;
      if (selected==-1) {
        if (sqrt(pow(mouseX-x, 2)+pow(mouseY-y, 2))<25) {
          selected = myId;
          return;
        }
      }
    }

    if (!UPDATES_ON_DRAW&&!drawsLines()) {
      return;
    }


    //Update position:
    float[] dispAdd = dispAddMemorySaver; //Sigma holder
    dispAdd[0] = 0;
    dispAdd[1] = 0;
    acc[0] = 0; 
    acc[1] = 0; //Acceleration calculation follows:
    for (int k = 0; k < distances.length; k++) {
      if (k==myId) {
        continue;
      }
      if (!nodes[k].exists) {
        continue;
      }
      float[] comp = componentsTo(nodes[k].pos);
      float compLength = sqrt(comp[0]*comp[0]+comp[1]*comp[1]);
      float dispar;
      float proximity = (compLength*(ZOOM/640.f)); //pixel to pd
      if (distances[k]>CURRENT_CUTOFF) {
        if (true) continue;
        dispar = 0;
        //If close proximity, move away.
        /*
        if (proximity < CURRENT_CUTOFF/2){
         dispar = .5;
         }
         */
      } 
      else {
        //Method 1
        dispar = (distances[k]-proximity) / (distances[k]+1);
        TOTAL_VARIANCE += dispar*dispar; //Unweighted.

        //Draw the line:
        drawLineTo(k, dispar);
      }
      //Weight for distance:
      //dispar /= (1+distances[k]);      
      //dispar = (exp(dispar-distances[k])-1.)/640.;

      if (proximity < 1e-5f) { //divide by zero catch
        dispAdd[0] = -dispar*1e-2f;
        dispAdd[1] = -dispar*1e-2f;
        //dispAdd[0] = dispar*1e-2;
        //dispAdd[1] = dispar*1e-2;
      } 
      else {
        dispAdd[0] = -dispar*(comp[0])/proximity;
        dispAdd[1] = -dispar*(comp[1])/proximity;
      }
      addVec(acc, dispAdd);
    }
    acc[0] /= mass;
    acc[1] /= mass;
    acc[0]=-vel[0]/timeStep*fluidFriction+acc[0];
    acc[1]=-vel[1]/timeStep*fluidFriction+acc[1];

    //Ok, move the particle, only if we're updating:

    if (UPDATES_ON_DRAW) {

      addVec(vel, acc[0]*timeStep, acc[1]*timeStep); 
      addVec(pos, vel[0]*timeStep, vel[1]*timeStep);
      float mag = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
      mag = min(min(mag, max(640, width/2)), max(640, height/2))/mag; //Bounding circle
      pos[0] *= mag;
      pos[1] *= mag;
    }
  }
  public float[] componentsTo(float[] other) {
    componentsMemSaver[0] = other[0]-pos[0];
    componentsMemSaver[1] = other[1]-pos[1];
    return componentsMemSaver;
  }
}
float[] dispAddMemorySaver = new float[2];
float[] componentsMemSaver = new float[2];
public void addVec(float[] a, float[] b) {
  a[0]+=b[0];
  a[1]+=b[1];
}
public void addVec(float[] a, float b1, float b2) {
  a[0]+=b1;
  a[1]+=b2;
}




public void TestFasta(int mode){

  FileDialog fd = new FileDialog(frame, "Please select your input file", FileDialog.LOAD);
  fd.show();
  File got = null;
  String toRead = "";
  if ( fd.getDirectory() != null && fd.getFile() != null) 
  {
    try {
      got = new File( fd.getDirectory() + File.separator + fd.getFile());
      if (!got.exists()){
        JOptionPane.showMessageDialog(frame, "You entered a nonexistant file:\n"+got);
        TestFasta(mode);
        return;
      }
      toRead = got.getAbsoluteFile().toString();
    } 
    catch (Throwable e){
      e.printStackTrace();
    };
    //println("You chose to open this file: " + chooser.getSelectedFile().getName());
  }  
  else {
    System.exit(0);
  }

  // String gotString = JOptionPane.ShowInputDialog("Please input a meaningful description of this graph.");
  String[] lines = loadStrings(toRead);
  File outtt = new File(got.toString()+".out");
  try {
    outtt.createNewFile();
    System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(outtt))));
  } 
  catch (Throwable e){
    e.printStackTrace();
    die("");
  };
  try {
    if (mode==0){
      Fastize(lines); 
      return;
    }
    if (mode==1){
      unFasta(lines,false); 
      return;
    }
    if (mode==2){
      removeDuplicates(lines); 
      return;
    }
    if (mode==3){
      unFasta(lines,true); 
      return;
    }
    if (mode==4){
      clustalwAnalyze(lines); 
      return;
    }
    if (mode==5){

      pdScore2Mat(lines);
      return;
    }

    if (mode==6){
      pdScoreMat(lines);
      return;
    }

    if (mode==7){
      seqSimScore(lines);
      return;
    }

    if (mode==8){
      pdScoreMultisearch(lines);
      return;
    }
  } 
  finally {
    System.out.flush();
    JOptionPane.showMessageDialog(frame,"O.K. The output has been stored to "+outtt);
  }
}
public void removeDuplicates(String[] lines){
  TreeSet ts = new TreeSet();
  for(int k = 0; k < lines.length; k++){
    ts.add(lines[k]);
  }
  Iterator i = ts.iterator();
  while(i.hasNext()){
    System.out.println(i.next());
  }
}

public String[][] unFastaSuper(String[] lines){
  ArrayList unFasta = new ArrayList();
  ArrayList headers = new ArrayList();
  StringBuilder building = null;
  for(int k = 0; k < lines.length+1; k++){
    boolean extraLine = k==lines.length;
    if (extraLine || lines[k].startsWith(">")){
      if (building!=null){
        String buildString = building.toString();
        if (buildString.matches("\\s+")){
          //Empty fasta sequence?
          throw new Error();
        } else {
          unFasta.add(buildString);
        }
      }
      building = new StringBuilder();
      if (extraLine) continue;
      headers.add(lines[k].substring(1));
    } 
    else {
      if (building==null){
        throw new Error();
      }
      building.append(lines[k]);
    }
  } 
  String[][] toRet = new String[2][unFasta.size()];
  for(int k = 0; k < unFasta.size(); k++){
    if (k < headers.size()){
      toRet[0][k]=(String)headers.get(k);
    } 
    else {
      toRet[0][k]="???";
    }
    toRet[1][k]=(String)unFasta.get(k);
    //System.out.println(k+" "+toRet[0][k]+" "+toRet[1][k]);
  }
  return toRet;
}
public void unFasta(String[] lines, boolean printHeaders){
  ArrayList unFasta = new ArrayList();
  ArrayList headers = new ArrayList();
  unFasta0(unFasta,headers,lines);
  for(int k = 0; k < unFasta.size(); k++){
    if (printHeaders){
      System.out.println(headers.get(k));
    } 
    else {
      System.out.println(unFasta.get(k));
    }
  }
}
public void unFasta0(ArrayList unFasta, ArrayList headers, String[] lines){
  String building = null;
  for(int k = 0; k < lines.length+1; k++){
    boolean extraLine = k==lines.length;
    if (extraLine || lines[k].startsWith(">")){
      if (building!=null) {
		  if (building.matches("\\s+")) {
	          throw new Error();
		  }
	      unFasta.add(building);
      }
      if (extraLine) continue;
      building = "";
      headers.add(lines[k].substring(1));
    } 
    else {
      building += lines[k];
    }
  } 
}
public void Fastize(String[] lines){
  int shortFormatThresh = 15;
  //Check if we can do sequence-headers
  boolean shortFormat = true;
  for(int k = 0; k < lines.length; k++){
    if (lines[k].length()>shortFormatThresh){
      shortFormat = false;
      break;
    }
  }
  for(int k = 0; k < lines.length; k++){
    if (lines[k].matches("\\s*")) continue;
    System.out.println(">"+(shortFormat?(""):(k+"-"))+lines[k].substring(0,min(lines[k].length(),shortFormatThresh)));
    for(int c = 0; c < lines[k].length(); c+=80){
      System.out.println(lines[k].substring(c,min(lines[k].length(),c+80)));
    }
  }
}
public void clustalwAnalyze(String[] lines){
  ArrayList validLines = new ArrayList(30000);
  int[] vals = new int[3];
  boolean inAlignment = false;
  for(int k = 0; k < lines.length; k++){
    if (lines[k].contains("Start of Pairwise alignments")){
      inAlignment = true;
    }
    if (lines[k].contains("Start of Multiple Alignment")){
      inAlignment = false;
    }
    if (!inAlignment){
      continue;
    }
    String[] doLine = lines[k].split("[^\\d]+");
    int got = 0;
    for(int p = 0; p < doLine.length; p++){
      if (!doLine[p].matches("\\s*")){
        vals[got]=new Integer(doLine[p]);
        got++;
        if (got==3){
          got = 0;
          validLines.add(vals);
          vals = new int[3];
        }
      }
    }
  }
  int Ncomp = validLines.size();
  int s = PApplet.parseInt((1+sqrt(1+8*Ncomp))/2); //AHAHAHHA fancy
  float[][] toFill = new float[s][s];
  for(int x = 0; x < s; x++){
    toFill[x][x]=0;
  }
  for(int k = 0; k < Ncomp; k++){
    vals = (int[])validLines.get(k);
    float value = min(20,1/pow(vals[2]/40f,2));
    toFill[vals[0]-1][vals[1]-1] = value;
    toFill[vals[1]-1][vals[0]-1] = value;
  }
  for(int y = 0; y < s; y++){
    for(int x = 0; x < s; x++){
      System.out.print(toFill[y][x]+" ");
    }
    System.out.println();
  }
}
public void saveNodeStates(String outputFile){
  try {
    File got = null;
    if (outputFile==null){
      //Save to file:
      FileDialog fd = new FileDialog(frame, "Please select where you would like to save the output (Include a file extension, example: .txt!)", FileDialog.SAVE);
      fd.show();
      String toRead = "";
      if ( fd.getDirectory() != null && fd.getFile() != null) 
      {
        try {
          got = new File( fd.getDirectory() + File.separator + fd.getFile());
          if (got.exists()){
            int gotOption = JOptionPane.showConfirmDialog(frame, "You entered an existant file:\n"+got+"\nOverwrite it?","Existing file.",JOptionPane.YES_NO_OPTION);
            if (gotOption == JOptionPane.YES_OPTION){
              //Just go on.
            } 
            else { 
              saveNodeStates(outputFile); //Null.
              return;
            }
          }
        } 
        catch (Throwable e){
          e.printStackTrace();
        };
        //println("You chose to open this file: " + chooser.getSelectedFile().getName());
      }  
      else {
        return; //Do nothing.
      }
    } 
    else {
      got = new File(outputFile);
    }

    toRead = got.getAbsoluteFile().toString();

    float PIX_TO_PD_FACTOR = (ZOOM/640.f);

    try {
      PrintWriter out = new PrintWriter(new FileWriter(got));
      for(int k = 0; k < nodes.length; k++){
        if (!nodes[k].exists) continue;
        out.println(">"+nodes[k].FastaLabel+" #COORDS:"+nodes[k].pos[0]*PIX_TO_PD_FACTOR+","+nodes[k].pos[1]*PIX_TO_PD_FACTOR);
        String seq = nodes[k].FastaSequence;
        for(int c = 0; c < seq.length(); c+=80){
          out.println(seq.substring(c,min(seq.length(),c+80)));
        }
      }
      out.close();

      out = new PrintWriter(new FileWriter(new File(toRead+".coords")));
      for(int k = 0; k < nodes.length; k++){
        if (!nodes[k].exists) continue;
        out.printf("%4d %-12s %8.3f %8.3f",k+1,nodes[k].FastaLabel,nodes[k].pos[0]*PIX_TO_PD_FACTOR,nodes[k].pos[1]*PIX_TO_PD_FACTOR);
        out.println();
      }
      out.close();
    } 
    catch (Throwable e){
      JOptionPane.showMessageDialog(frame,"I could not write to the file \n"+got+"\n due to "+e.toString(),"Write Fasta Error",JOptionPane.ERROR_MESSAGE); 
    }
  } 
  finally { 
    //Inscript check
    if (!iS){
      g.endDraw();
      frame.setSize(new Dimension(frame.getWidth()-10,frame.getHeight()-10));
      try {
        Thread.sleep(100);
      } 
      catch (Throwable e){
      };
      frame.setSize(new Dimension(frame.getWidth()+10,frame.getHeight()+10));
      g.beginDraw();
    }
  }
}

public void GetCustomDistances(String customLoc){
  try {
    /*
    Scanner input = new Scanner(new BufferedInputStream(openStream(customLoc),262144));
     int lineNum = 0;    
     while(input.hasNextLine()){
     input.nextLine(); 
     lineNum++;
     }
     input.close();
     System.out.println("Second pass through custom matrix...");
     input = new Scanner(new BufferedInputStream(openStream(customLoc),262144));
     distanceMatrix = new float[lineNum][lineNum];
     for(int k = 0; k < lineNum; k++){
     //String[] thisLine = input.nextLine().split("\\s+");
     for(int p = 0; p < lineNum; p++){
     distanceMatrix[k][p]=new Float(input.next());
     }
     input.nextLine();
     }
     input.close();
     */
    byte[] buffer = new byte[0x40000]; //256kb at a time
    for(int actuallyFillBuffer = 0; actuallyFillBuffer < 2; actuallyFillBuffer++){
      //Used to be called openStream...
      InputStream is = openStream(customLoc);
      try { 
        int read = -1;
        int lineNum = 0;
        boolean fillBuffer = actuallyFillBuffer==1;
        int linex = 0;
        boolean stillInLineReturn = false;
        String leftOver = "";
        while((read=is.read(buffer))!=-1){
          String k = leftOver+new String(buffer,0,read);
          int kLenNoTail = -1;
          int y = 0;
          for(y = k.length()-1; y >=0; y--){
            if(Character.isWhitespace(k.charAt(y))){
              kLenNoTail = y+1;
              break;
            }
          }
          for(y = 0; y < kLenNoTail; y++){
            if (isLineEndChar(k.charAt(y))){
              if (stillInLineReturn){
                continue;
              }
              stillInLineReturn = true;
              if (fillBuffer&&lineNum>0&&linex < distanceMatrix.length){
                throw new RuntimeException("Assertion error: not square matrix.");
              }
              linex = 0;
              lineNum++;
              while ((y+1<kLenNoTail)&&isLineEndChar(k.charAt(y+1))){
                y++;             
              }
            } 
            else {
              stillInLineReturn = false;
              if (fillBuffer){
                //Parse float.
                if (!Character.isWhitespace(k.charAt(y))){
                  int begin = y;
                  int end = y;
                  while ((y+1<kLenNoTail)&&!Character.isWhitespace(k.charAt(y+1))){
                    y++;        
                    end = y;
                  }
                  distanceMatrix[lineNum][linex++] = new Float(k.substring(begin,end+1));
                  /*
                  if (distanceMatrix[lineNum][linex-1]==0){
                    //System.out.println("Zero"+" "+begin+" "+(end+1));
                    if (lineNum!=linex-1){
                      System.out.println("Bad!"+" "+lineNum+" "+(linex-1));
                    }
                  }
                  */
/*      
                  if (linex == 1){
                    System.out.println(lineNum+" "+distanceMatrix[lineNum][0]);
                  }
*/              
                }
              }  
            }
          }
          leftOver = k.substring(kLenNoTail,k.length());
        }
        if (!fillBuffer){
          System.out.println(lineNum);
          distanceMatrix = new float[lineNum][lineNum];
          lineNum = 0;
        }
      }
      finally {
        is.close();
      }
    }
  } 
  catch (Throwable e){
    e.printStackTrace();
  }
}
private static boolean isLineEndChar(char c){
  if (c > '\r'){
    return false;
  }
  return c=='\n'||c=='\r';
}
public void pdScore2Mat(String[] lines){
  int wSize = -1;
  String input =  JOptionPane.showInputDialog(frame,"Input the Window Size(integer):");
  try {
    wSize = new Integer(input);
  } 
  catch (Throwable e){
    JOptionPane.showMessageDialog(frame,("Non integer input."),"You entered a non-integer window size.",JOptionPane.ERROR_MESSAGE);
    return;
  }
  for(int k = 0; k < lines.length; k++){
    for(int p = 0; p < lines.length; p++){
      float val = PDScore2(lines[p],lines[k],wSize);
      System.out.print(val+" ");
    }
    System.out.println();
  }
}

public void pdScoreMat(String[] lines){
  for(int k = 0; k < lines.length; k++){
    for(int p = 0; p < lines.length; p++){
      float val = PDScore(lines[p],lines[k]);
      System.out.printf("%.3f ",val);
    }
    System.out.println();
    System.err.println(k+" out of "+lines.length);
  }
}

public void seqSimScore(String[] lines){
  ArrayList unFasta = new ArrayList();
  ArrayList headers = new ArrayList();
  unFasta0(unFasta,headers,lines);

  for(int k = 0; k < unFasta.size(); k++){
    for(int p = 0; p < unFasta.size(); p++){
      System.out.printf("%.3f ",seqSimDist((String)unFasta.get(k),(String)unFasta.get(p)));   
    }
    System.out.println();
    System.out.flush();
  }
}
public float seqSimDist(String a, String b){
  float diff = 0;
  for(int k = 0; k < a.length(); k++){
    if (k >= b.length() || a.charAt(k)!=b.charAt(k)){
      diff++;
    }
  }
  return diff / sqrt(a.length()*a.length()+b.length()*b.length()) * 25;
}

public void pdScoreMultisearch(String[] lines){
  FileDialog fd = new FileDialog(frame, "Choose the short list of peptides (to match against)", FileDialog.LOAD);
  fd.show();
  File got = null;
  String toRead = "";
  if ( fd.getDirectory() != null && fd.getFile() != null) 
  {
    try {
      got = new File( fd.getDirectory() + File.separator + fd.getFile());
      if (!got.exists()){
        JOptionPane.showMessageDialog(frame, "You entered a nonexistant file:\n"+got);
        return;
      }
      toRead = got.getAbsoluteFile().toString();
    } 
    catch (Throwable e){
      e.printStackTrace();
    };
    //println("You chose to open this file: " + chooser.getSelectedFile().getName());
  }  
  else {
    System.exit(0);
  }  
  String[] shortlist = loadStrings(toRead);

  double thresh = -1;
  String input =  JOptionPane.showInputDialog(frame,"Input the PD similarity required:");
  try {
    thresh = new Double(input);
  } 
  catch (Throwable e){
    JOptionPane.showMessageDialog(frame,("Non numeric input."),"You entered a non-numeric threshold.",JOptionPane.ERROR_MESSAGE);
    return;
  }

  //May print duplicates
  for(int k = 0; k < lines.length; k++){
    String peptide = lines[k];
    for(int p = 0; p < shortlist.length; p++){
      if (PDScore(peptide,shortlist[p])<thresh){
        System.err.println(peptide+" "+shortlist[p]+" "+thresh);
        System.out.println(peptide);
        break;
      }
    }
  }

}

public float PDScore(String a, String b){
  a = a.trim().toUpperCase();
  b = b.trim().toUpperCase();
  //PRE: a shorter than B.
  String smaller = a.length()<b.length()?a:b;
  String bigger = a.length()<b.length()?b:a;
  a = smaller;
  b = bigger;
  //Done.
  float toMin = Float.MAX_VALUE;
  int k = 0;
  for(k = 0; k <= bigger.length()-smaller.length(); k++){
    toMin = min(toMin, PDScore0(a,b.substring(k,k+a.length())));
  }
  //System.out.println(a+" "+b+" "+(k-1)+" "+toMin);
  return toMin;
}
public float PDScore0(String a, String b){
  float totSum = 0;
  float partSum = 0;
  float pval1 = 0;
  float pval2 = 0;
  int att = 0;
  int k = 0;
  for(k = 0; k < b.length(); k++){
    partSum = 0;
    for(att = 0; att < 5; att++){
      pval1 = pdVal(a.charAt(k),att);
      pval2 = pdVal(b.charAt(k),att);
      pval1 = pval1-pval2;
      partSum = partSum + pval1*pval1;
    }
    totSum += sqrt(partSum);
  }
  return totSum/a.length();
}
char[] list1 = new char[]{
  'A',	'R',	'N',	'D',	'C',	'Q',	'E',	'G',	'H',	
  'I',	'L',	'K',	'M',	'F',	'P',	'S',	'T',	'W',	'Y',	'V',	'-'};
int[] list;
{
  list = new int[128];
  for(int k = 0; k < list1.length; k++){
    list[list1[k]]=k;
  }
}
float[][] atts = new float[][]{
  {
    0.354f,	7.573f,	11.294f,	13.42f,	-5.846f,	6.599f,	9.788f,	9.655f,	1.019f,	-15.634f,	-11.825f,	10.762f,	-10.585f, -14.571f, 7.662f, 8.813f,	3.012f,	-13.11f,	-6.245f,	-12.135f, 0,  }
  ,
  {
    3.762f,	-10.135f,	1.067f,	-1.6f,	4.885f,	-5.166f,	-7.861f,	15.778f,	-4.969f,	1.993f,	0.505f,	-9.517f,	-3.959f,	-0.646f,	8.029f,	6.682f,	4.127f,	-5.222f,	-1.6f,	3.818f,	0,  }
  ,
  {
    -11.036f,	2.486f,	2.718f,	-0.325f,	1.626f,	-0.697f,	-7.318f,	-0.558f,	0.953f,	-2.045f,	-6.157f,	-1.022f,	-3.601f,	1.673f,	9.456f,	-0.348f,	-0.348f,	9.038f,	9.874f,	-4.345f,	0,  }
  ,
  {
    -0.649f,	-4.291f,	1.963f,	3.742f,	9.397f,	0.582f,	2.611f,	0.299f,	4.657f,	-3.243f,	-4.557f,	-5.405f,	5.339f,	-0.033f,	-3.576f,	-1.131f,	-2.195f,	1.38f,	-1.597f,	-3.26f,	0,  }
  ,
  {
    2.828f,	-5.687f,	-0.859f,	2.437f,	-5.843f,	-1.75f,	4.734f,	1.656f,	-0.328f,	-1.672f,	3.219f,	-0.422f,	1.203f,	3.25f,	6,	-3.062f,	-4.281f,	4.64f,	-1.422f,	-4.672f,	0,  }
};

public float pdVal(char one, int att){
  int whichChar = list[(int)one];
  /*
   for(int k = 0; k < list.length; k++){
   if (one == list[k]){
   whichChar = k;
   }
   }
   */
  return atts[att][whichChar];
}

/**
 * Computes the PD2 score between sequences A and B. Both may be of variable lengths.
 **/
public float PDScore2(String a, String b, int wSize){
  if (a.length() > b.length()){
    String tmp = a;
    a = b;
    b = tmp;
  }

  boolean[] bFlagged = new boolean[b.length()];
  //the last wSize-1 of bFlagged can never possibly be truth'ed.
  int matches = 0;
  for(int w = 0; w < a.length()-wSize; w++){
    //Check if b contains this string:
    String toFind = a.substring(w,w+wSize);
    for(int search = 0; search < b.length(); search++){
      //TODO: search for approximate matches 
      int value = b.indexOf(toFind,search);
      if (value==-1){
        break;
      }
      if (!bFlagged[value]){
        bFlagged[value] = true;
        matches++;
        break;
      }
      search = value+1;
    }
  }
  return 20* sqrt( 1 - (matches / (float) (a.length()-wSize)));
}
public void getCoordsFile(float[][] coords, String file){
  String[] strings = loadStrings(file);
  for(int k = 0; k < strings.length; k++){
    String[] coord = strings[k].split("\\s+");
    if(coord.length < 5){
      continue;
    }
    coords[k][0] = new Float(coord[3]);
    coords[k][1] = new Float(coord[4]);
    float PIX_TO_PD_FACTOR = (ZOOM/640.f);
    coords[k][0]/=PIX_TO_PD_FACTOR;
    coords[k][1]/=PIX_TO_PD_FACTOR;
  }
}
/**
An algorithm that rotates / flips an input set (toModify) of 2D points to best match another set (best), using tmpMem as a required midway buffer. 
**/
public void rotateToBest(float[][] best, float[][] toModify, float[][] tmpMem, String outFile, PrintWriter dbg){
  float bestMag = Float.MAX_VALUE;
  int needsAFlip = 0;
  float andTheta = 0;
  for(int doIt = 0; doIt < 2; doIt++){
    float theta = calcBestRotation(best,toModify);
    for(float offTheta = 0; offTheta < 1; offTheta++){ //nO ADD PIE!@
      rotateInto(toModify,tmpMem,theta);
      float mag = sumDeviant(best,tmpMem);
      if (mag < bestMag){
        bestMag = mag;
        andTheta = theta;
        needsAFlip = doIt;
      }
      theta+=PI;
    }
    //For second time, try flipping the x-axis of tomodify.
    for(int k = 0; k < toModify.length; k++){
      toModify[k][0]=-toModify[k][0];
    }
  }
  //if the flip was better, flip to get there..
  if (needsAFlip==1){
    for(int k = 0; k < toModify.length; k++){
      toModify[k][0]=-toModify[k][0];
    }
  }
  rotateInto(toModify,tmpMem,andTheta);
  dbg.printf(" %8.3f %-5d %8.3f",andTheta,needsAFlip,bestMag);
  dbg.println();
  //Ok. tmpmem holds the ideal rotation, write it out:
  for(int k = 0; k < nodes.length; k++){
    nodes[k].pos = tmpMem[k];
  }
  saveNodeStates(outFile);
}
public void rotateInto(float[][] one, float[][] dst, float theta){
  float ct = cos(theta); 
  float st = sin(theta);
  for(int k = 0; k < one.length; k++){
    dst[k][0] = ct*one[k][0]+st*one[k][1];
    dst[k][1] = -st*one[k][0]+ct*one[k][1];
  }
}
/**
Calculates the best rotation to fit one onto two. Note that each point of one is calculated as being the 'partner' of the corresponding indexed point in two. 
A more complex algorithm could be developed that rotates generalized point cloud A onto generalized point cloud B, by somehow choosing the partners on the fly. Here, we know which points
correspond (they are data points).
**/
public float calcBestRotation(float[][] one, float[][] two){
  float[] q = new float[4];
  for(int k = 0; k < one.length; k++){
    q[0]+=one[k][0]*two[k][0];
    q[1]+=one[k][0]*two[k][1];
    q[2]+=one[k][1]*two[k][0];
    q[3]+=one[k][1]*two[k][1];
  }
  float arcus = atan2((q[1]-q[2]),(q[0]+q[3]));
  return arcus;
}
public float sumDeviant(float[][] one, float[][] two){
  float toRet = 0;
  float dx, dy;
  for(int k = 0; k < one.length; k++){
    dx = one[k][0]-two[k][0];
    dy = one[k][1]-two[k][1];
    toRet += sqrt(dx*dx+dy*dy);
  }
  return toRet;
}
  static public void main(String[] passedArgs) {
    //If this returns the right value, then HIDPI scaling will happen if the
    //system default LAF is the GTK one.
    //System.out.println(Toolkit.getDefaultToolkit().getScreenResolution());
    String[] appletArgs = new String[] { "ClusterSimulation" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
