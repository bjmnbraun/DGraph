package com.dgraph;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.GridLayout;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.TreeSet;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextArea;

import com.pdscore.PDScore;

import processing.core.PApplet;
import processing.core.PFont;

@SuppressWarnings("serial")
public class DGraph extends PApplet {
	String removed = "";

	boolean PRESETUP = false;
	int RUN_MODE = -1;
	int INPUT_MODE = -1;
	String toRead;
	String customLoc;
	boolean noninteractive = false;
	boolean singleFramePdf = false;

	public void start() {
		//Test here if we have write priviledges
		String test = ".test";
		try {
			FileOutputStream fos = new FileOutputStream(test);
			fos.write(2);
			fos.close();
		} 
		catch (Throwable e){
			JOptionPane.showMessageDialog(frame,"Please run from a directory that has file-write access.\n","Write-Protected Directory",JOptionPane.ERROR_MESSAGE);
			exit();
		}

		doCommandLine();
		if (singleFramePdf){
			PRESETUP = true;
		}
		if (!PRESETUP){
			PRESETUP = true;
			//Show a dialog with options for RUN_MODE
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
			holder.setVisible(true);
			boolean shouldRun = JOptionPane.showOptionDialog(
					holder,
					pane,
					"Welcome! What would you like to run?",
					JOptionPane.OK_CANCEL_OPTION,
					JOptionPane.PLAIN_MESSAGE,null,null,null) == JOptionPane.OK_OPTION;
			if (!shouldRun) {
				exit();
			}
			holder.setVisible(false);

			int choice = 0;
			for(int k = 0; k < options.length; k++){
				if (options[k].isSelected()){
					choice = k;
					break;
				}
			}
			RUN_MODE = choice-1; //0=-1
		}

		if (RUN_MODE != -1){
			TestFasta(RUN_MODE);
			exit();
		}
		
		destroySimulation();
		setupSimulation();
		
		//Now do PApplet's start to get things rolling.
		super.start();
	}

	//This is actually run on the first "frame." A call to size triggers it to be called _again_ (!)
	public void setup() {
		size(800,800,P2D);
		if (!frame.isResizable()){
			frame.setResizable(true);
		}
		txt = createFont("Vera.ttf",90,true);
		frameRate(60);
	}

	long lastRenderTime = 0;
	public void draw(){
		if (!frame.isResizable()){
			frame.setResizable(true);
		}
		//TODO do this in a more sensible way...
		//drawSimulation(true);
		if (keyPressed && key=='p' && (System.nanoTime() - lastRenderTime) > 0.5e9){
			renderThisFrame = true;
			lastRenderTime = System.nanoTime();
		}
		drawOrRecordSimulation();
	}

	float[] shadings;
	int[] nodeConnections;
	String renderUrl = null;
	//Output directory for drawShaded
	String renderDir = null;
	public void setupShaded() {
		shadings = new float[nodes.length]; //defaults to 0
		nodeConnections = new int[nodes.length];

		String toRead = null;
		String[] commands = null;
		if (singleFramePdf){
			commands = new String[0];
		}
		if (commands==null){
			FileDialog fd = new FileDialog(frame, "Open the coloring file", FileDialog.LOAD);
			//while ( fd.getDirectory() == null && fd.getFile() == null) 
			//{
			fd.show();
			//}
			try {
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

		for(int k = 0; k < nodes.length; k++){
			if (!nodes[k].display) continue;
			nodes[k].shadeColor = lerpColor(color(0,0,0),color(255,0,0),constrain(shadings[k],0,1)); //erased after 1 use, by the way.
			nodes[k].nodeConnection = nodeConnections[k];
		}
	}
	boolean renderThisFrame = false;
	public void drawOrRecordSimulation(){
		String renderLocation = null;
		if (renderThisFrame){
			renderLocation = renderUrl!=null?renderUrl:((renderDir!=null?renderDir+"/":"")+frameCount);
			beginRecord(PDF, renderLocation + ".pdf");
		}
		drawSimulation(!renderThisFrame);
		if (renderThisFrame){
			endRecord();
			//Try to convert pdf to png with imagemagick, if installed:
			ProcessBuilder pb = new ProcessBuilder("magick","-density","300",renderLocation + ".pdf",/*"-resize","25%",*/renderLocation + ".png");   
			pb.redirectErrorStream(true);   
			try {   
				Process p = pb.start();   
				BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));   
				String line = null;   
				while((line=br.readLine())!=null){   
					System.out.println(line);   
				}   
				int returnCode = p.waitFor();
			} catch(Throwable e) {
				e.printStackTrace();
			} 
			renderThisFrame = false;
			if (singleFramePdf){
				exit();
			}
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
			exit();
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
		noninteractive = true; //inScript flag.

		File outputDir = new File(outputDirectory).getAbsoluteFile();
		if (!outputDir.isDirectory() || !outputDir.exists()){
			batchError("-Dout must be a pre-existing directory."+outputDir);
		}
		File inputFasta = new File(inFasta).getAbsoluteFile();
		if (!inputFasta.exists()){
			batchError("-Din referred to bad/no file\n"+inFasta);
		}
		toRead = inFasta; //Read by clusterSimulation.
		INPUT_MODE = 0;
		if (customMatrix!=null){
			File inputCustomMatrix = new File(customMatrix);
			if (!inputCustomMatrix.exists()){
				batchError("-DinCustom referred to bad/no file\n"+inputCustomMatrix);
			}
			customLoc = customMatrix;
			INPUT_MODE = 1;
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
			renderUrl = outputDir.toString()+File.separator+filePrefix+".pdf";
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
			float[][] best = new float[numNodes][2];
			getCoordsFile(best, filesWrote[bestRunIndex]+".fasta.coords");
			float[][] otherRun = new float[numNodes][2];
			float[][] out = new float[numNodes][2];
			ana.println("Start of rotations:");
			for(int run = 0; run < numTimes; run++){
				if (run==bestRunIndex) continue; //TODO can probably remove this
				getCoordsFile(otherRun, filesWrote[run]+".fasta.coords");
				//So, modify filesWrote[run] so that it's better.
				ana.printf("%10d",run);
				PointCloudUtil.transformToMatch2D(best, otherRun, out,  ana);
				//Ok. out holds the RMSD minimized transformed otherRun, write it out:
				for(int k = 0; k < nodes.length; k++){
					nodes[k].pos = out[k];
				}
				saveNodeStates(filesWrote[run]+".fasta");
			}
			ana.close();
		} 
		catch (Throwable e){
			e.printStackTrace();
			batchError("Error in analysis:"+e.toString());
		}
		exit();
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
		} 
		return "data.txt";
	}
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
		boolean useColoringFile = false;
		if (INPUT_MODE == -1) {
			String [] optionsStr = { 
					"Fasta sequences, generate PD scores",
					"Fasta sequences, use custom distance matrix",
					"Just labels (one per line), use custom distance matrix",
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
			JCheckBox useColoringFileCB = new JCheckBox("Use coloring file?");
			pane.add(useColoringFileCB);
			options[0].setSelected(true);
			useColoringFileCB.setSelected(false);

			Frame holder = new Frame("DGraph - Initializing.");
			holder.setVisible(true);
			boolean shouldRun = JOptionPane.showOptionDialog(
					holder,
					pane,
					"Select input data type",
					JOptionPane.OK_CANCEL_OPTION,
					JOptionPane.PLAIN_MESSAGE,null,null,null) == JOptionPane.OK_OPTION;
			if (!shouldRun) {
				exit();
			}
			holder.setVisible(false);

			INPUT_MODE = 0;
			for(int k = 0; k < options.length; k++){
				if (options[k].isSelected()){
					INPUT_MODE = k;
					break;
				}
			}
			useColoringFile = useColoringFileCB.isSelected();
		}

		if (INPUT_MODE == 0) {
			//Fasta, generate PD scores
			if (toRead == null) {
				toRead = pickUserFile("Open a Fasta File");
			}
			String[][] fastaInfo = loadFasta(toRead);
			initNodesFromFasta(fastaInfo);
			initDistancesWithPDScore(fastaInfo);
		}
		if (INPUT_MODE == 1) {
			//Fasta, generate 
			if (toRead == null) {
				toRead = pickUserFile("Open a Fasta File");
			}
			String[][] fastaInfo = loadFasta(toRead);
			initNodesFromFasta(fastaInfo);
			if (customLoc == null) {
				scoreName = JOptionPane.showInputDialog(frame, "Enter the name of your custom score function");
				customLoc = pickUserFile("Select " + scoreName + " matrix file");
			}
			loadDistancesFromFile(customLoc);
		}
		//Make a directory in the same directory as toRead.
		renderDir = toRead + "-DGraph";
		new File(renderDir).mkdir();
		if (INPUT_MODE == 2) {
			//Custom labels not yet implemented
			throw new RuntimeException("Not yet implemented.");
		}
		if (useColoringFile) {
			setupShaded();
		}
		removeIsolatedNodes();
		//TODO should be a single button to print PDF + PNG.

		//if (!iS) g.beginDraw();
	}

	private void loadDistancesFromFile(String customLoc) {
		try {
			byte[] buffer = new byte[0x40000]; //256kb at a time
			for(int actuallyFillBuffer = 0; actuallyFillBuffer < 2; actuallyFillBuffer++){
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
								if (fillBuffer&&lineNum>0&&linex < nodes.length){
									throw new RuntimeException("Row " + lineNum + " did not match # of nodes (" + nodes.length+")");
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
										if (nodes[lineNum].distances == null) {
											nodes[lineNum].distances = new float[nodes.length];
										}
										nodes[lineNum].distances[linex++] = new Float(k.substring(begin,end+1));
									}
								}  
							}
						}
						leftOver = k.substring(kLenNoTail,k.length());
					}
					if (!fillBuffer){
						System.out.println(lineNum);
						lineNum = 0;
					}
				}
				finally {
					is.close();
				}
			}
		} 
		catch (Throwable e){
			JOptionPane.showMessageDialog(frame, "Error reading distance matrix: " + e.toString(), "Error!", JOptionPane.ERROR_MESSAGE);
			exit();
		}
	}
	private void initDistancesWithPDScore(String[][] fastaInfo) {
		scoreName = "PD";

		long lastDraw = System.nanoTime();
		for (int k = 0; k < nodes.length; k++) {
			float[] dist = new float[nodes.length];
			for (int p = 0; p < nodes.length; p++) {
				if (p==k) continue;
				dist[p]=(float)PDScore.PDScore(fastaInfo[1][p], fastaInfo[1][k]);
			}
			nodes[k].distances = dist;
			nodes[k].myId = k;
			nodes[k].display = true;
			//System.out.println("Processed: "+lines[k]);
			if (!noninteractive) if ((System.nanoTime()-lastDraw)>.5e9f) {
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
	}
	private void removeIsolatedNodes() {
		float dontShowCutoff = REMOVAL_ISLAND_FINDER;
		boolean hasLessThanCutoff = false;
		for (int k = 0; k < nodes.length; k++) {
			float[] dist = nodes[k].distances;
			for (int p = 0; p < nodes.length; p++) {
				if (p==k) continue;
				if (dist[p] < dontShowCutoff /**&& dist[p]!=0**/)
					hasLessThanCutoff = true;
				if (!hasLessThanCutoff) {
					removed+=" ,"+(k+1);
					nodes[k].display = false;
				}
			}
		}
	}

	private void initNodesFromFasta(String[][] fastaInfo) {
		nodes = new Node[fastaInfo[0].length];
		for (int k = 0; k < nodes.length; k++) {
			nodes[k] = new Node();
		}
		removed = "";

		//if (!iS) noLoop();
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
		}


	}
	private String[][] loadFasta(String toRead) {
		String[] lines = loadStrings(toRead);
		String[][] fastaInfo = null;
		try {
			fastaInfo = unFastaSuper(lines);
		} catch (Throwable e) {
			JOptionPane.showMessageDialog(frame, "Error reading fasta file: " + e.toString(), "Error!", JOptionPane.ERROR_MESSAGE);
			exit();
		}
		return fastaInfo;
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
			if (!got.display) continue;
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
	private int NodeDotSize = 3; //rectangle of this radius at each node
	private String underMouse;
	int simulationTicks = 0;

	private boolean UPDATES_ON_DRAW = true; //Normally true, but set to false on the first frame in the GUI.
	private boolean wantsFrameBreak = false;
	private boolean simStarted = false;
	public void drawSimulation(boolean controlInfo) {
		background(255);
		UPDATES_ON_DRAW = simStarted; //GUI uses this so the user can just look at it for a second.
		underMouse = "";
		textFont(txt);
		textSize(txtFontSizeMedium * txtFontSizeScaling);
		noStroke();
		fill(0);
		if (controlInfo) {
			textAlign(LEFT, TOP);
			text("\"Islands\" Removed by "+REMOVAL_ISLAND_FINDER+" PD threshold: "+removed.substring(min(removed.length(), 2)), 2, controlInfo?20:5, width, 80);
			text("Noncomparable Seperation Factor:"+CURRENT_CUTOFF, 2, 4);
			text("Displaying:"+SHOW_LINES, 2, 100);
			textAlign(RIGHT, BOTTOM);
			if (!UPDATES_ON_DRAW) {
				fill(255, 0, 0);
				text("Simulation stopped. Press 'y' to start again.", width, 0);
			}
			fill(0);
			float ypos = height - textDescent();
			float ydelta = textAscent() + textDescent();
			float xpos  = width - 2;
			text("'h' to switch on/off lines", xpos, ypos -= ydelta);
			text("'c' to turn off labels", xpos, ypos -= ydelta);
			text("'s' to show sequence numbering", xpos, ypos -= ydelta);
			text("'a' to show sequences", xpos, ypos -= ydelta);
			text("'b' to show sequence headers", xpos, ypos -= ydelta);
			text("'9' zooms out", xpos, ypos -= ydelta);
			text("'0' zooms in", xpos, ypos -= ydelta);
			text("'m' saves the current state", xpos, ypos -= ydelta);
			text("'p' renders current state", xpos, ypos -= ydelta);
			text("'j' sets the physics options", xpos, ypos -= ydelta);
			text("'d' disable nodes", xpos, ypos -= ydelta);
			text("'w' increases dot-size", xpos, ypos -= ydelta);
			text("'q' decreases dot-size", xpos, ypos -= ydelta);
			text("'r' or 'e' re-run the simulation. 'e' only scatters the points, 'r' does more.", xpos, ypos -= ydelta);
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
					nodes[o].display = false;
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
			if (!nodes[k].display) continue;
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
		boolean display;
		int shadeColor = -1;
		int nodeConnection = -1;
		String FastaLabel;
		String FastaSequence;
		int myId;
		public void draw() {
			if (!display) {
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
			if (!display) {
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
			if (noninteractive) return;
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
			if (!display) {
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
				if (!nodes[k].display) {
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
			exit();
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
				ClustalwUtil.clustalwAnalyze(lines); 
				return;
			}
			if (mode==5){

				pdScore2Mat(lines);
				return;
			}

			if (mode==6){
				PDScore.pdScoreMat(lines);
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
				float val = (float) PDScore.PDScore2(lines[p],lines[k],wSize);
				System.out.print(val+" ");
			}
			System.out.println();
		}
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
			exit();
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
				if (PDScore.PDScore(peptide,shortlist[p])<thresh){
					System.err.println(peptide+" "+shortlist[p]+" "+thresh);
					System.out.println(peptide);
					break;
				}
			}
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
					if (!nodes[k].display) continue;
					out.println(">"+nodes[k].FastaLabel+" #COORDS:"+nodes[k].pos[0]*PIX_TO_PD_FACTOR+","+nodes[k].pos[1]*PIX_TO_PD_FACTOR);
					String seq = nodes[k].FastaSequence;
					for(int c = 0; c < seq.length(); c+=80){
						out.println(seq.substring(c,min(seq.length(),c+80)));
					}
				}
				out.close();

				out = new PrintWriter(new FileWriter(new File(toRead+".coords")));
				for(int k = 0; k < nodes.length; k++){
					if (!nodes[k].display) continue;
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
			if (!noninteractive){
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


	private static boolean isLineEndChar(char c){
		if (c > '\r'){
			return false;
		}
		return c=='\n'||c=='\r';
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

	static public void main(String[] passedArgs) {
		String[] appletArgs = new String[] { "com.dgraph.DGraph" };
		if (passedArgs != null) {
			PApplet.main(concat(appletArgs, passedArgs));
		} else {
			PApplet.main(appletArgs);
		}
	}
}
