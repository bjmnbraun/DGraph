package com.dgraph;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.Frame;
import java.awt.GridLayout;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.time.Instant;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;

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
	String inputFileName;
	String customLoc;
	boolean headless = false;
	int scriptRenderFrames = 0;
	boolean SHOW_HELP_TEXT = true;
	int VIEW_WIDTH = 800;
	int VIEW_HEIGHT = 800;

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
			System.exit(1);
		}

		doCommandLine();
		if (scriptRenderFrames > 0){
			PRESETUP = true;
		}
		while (!PRESETUP){
			PRESETUP = true;
			//Show a dialog with options for RUN_MODE
			String [] optionsStr = { 
					"Run DGraph", 
					"Utility: Sequence list -> Fasta", 
					"Utility: Fasta -> Sequence List", 
					"Utility: Remove Duplicates from Sequence List", 
					"Utility: Fasta -> List of headers",
					"Utility: Clustalw output -> Alignment score distance matrix. Note: Copy+paste whole clustalw output, especially 'pairwise alignment' section",
					"Utility: Fasta -> PD Score matrix: sliding window approach",
					"Utility: Fasta -> PD Score matrix: assume pre-aligned sequences",
					"Utility: DNA Fasta -> DNA Sequence Similarity matrix",
					"Utility: Find all pairs of aligned peptides in a list with PD Score under a threshold",
					"Utility: Jalview all pairwise sequeqnce alignments output -> Parwise alignment score list",
					"Utility: Reorder fasta sequences to match another file",
					"Utility: Fasta -> Line-numbered sequence headers",
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
			holder.setSize(200,100);
			holder.setLocationRelativeTo(null);
			holder.setVisible(true);
			boolean shouldRun = JOptionPane.showOptionDialog(
					holder,
					pane,
					"Select operation",
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
			RUN_MODE = choice-1; //0=-1
			
			if (RUN_MODE != -1){
				runUtility(RUN_MODE);
				//Go again.
				PRESETUP = false;
			}
		}

		headless = false;
		//Now do PApplet's start to get things rolling.
		super.start();
	}

	//This is actually run on the first "frame." A call to size triggers it to be called _again_ (!)
	public void setup() {
		//This throws an exception and causes us to reenter once "frame" exists
		size(VIEW_WIDTH,VIEW_HEIGHT,P2D);
		txt = createFont("Vera.ttf",90,true);
		frameRate(60);
		if (!frame.isResizable()){
			frame.setResizable(true);
		}
	}

	long lastRenderTime = 0;
	boolean needsSetupSimulation = true;
	public void draw(){
		if (!frame.isResizable()){
			frame.setResizable(true);
		}
		//Frame exists here.
		if (needsSetupSimulation) {
			destroySimulation();
			setupSimulation();
			needsSetupSimulation = false;
		}
		if (keyPressed && key=='s' && (System.nanoTime() - lastRenderTime) > 0.5e9){
			renderThisFrame = true;
			lastRenderTime = System.nanoTime();
		}
		drawOrRecordSimulation();
		if (!simStarted && !renderThisFrame) {
			noLoop();
		}	
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
		if (scriptRenderFrames > 0){
			commands = new String[0];
		}
		if (commands==null){
			FileDialog fd = new FileDialog(frame, "Open the coloring file", FileDialog.LOAD);
			//while ( fd.getDirectory() == null && fd.getFile() == null) 
			//{
			showFDRememberDirectory(fd);
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
			renderLocation = renderSaveLocationBase(); 
			if (scriptRenderFrames > 0) {
				if (scriptRenderFrames == 2) {
					showLabels = SHOW_LABELS_NAME;
					renderLocation += "-names";
				}
				if (scriptRenderFrames == 1) {
					showLabels = SHOW_LABELS_NUMBER;
					renderLocation += "-numbered";
				}
			}
			beginRecord(PDF, renderLocation + ".pdf");
			smooth();
		}
		drawSimulation(!renderThisFrame);
		if (renderThisFrame){
			endRecord();
			noSmooth();
			//Try to convert pdf to png with imagemagick, if installed:
			/*
			ProcessBuilder pb = new ProcessBuilder("magick","-density","300",renderLocation + ".pdf",renderLocation + ".png");   
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
			*/
			save(renderLocation + ".png");
			if (scriptRenderFrames > 0) {
				scriptRenderFrames--;
				if (scriptRenderFrames == 0) {
					System.exit(0);
				}
				//Keep going with renderThisFrame still true.
			} else {
				renderThisFrame = false;
			}
		}
	}
	private String renderSaveLocationBase() {
		//https://www.dotnetperls.com/filename-date-java
		// Get a Calendar and set it to the current time.
        Calendar cal = Calendar.getInstance();
        cal.setTime(Date.from(Instant.now()));

        // Create a filename from a format string.
        // ... Apply date formatting codes.
        String datestring = String.format("%1$tY-%1$tm-%1$td-%1$tk-%1$tM-%1$tS", cal);
		return renderUrl!=null?renderUrl:((renderDir!=null?renderDir+"/":"")+datestring);
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
			System.exit(0);
		}
		boolean doAutomation = false;
		String numRuns = System.getProperty("runs");
		String outputDir = System.getProperty("out");
		String filePref = System.getProperty("prefix");
		String numIters = System.getProperty("iters");
		String params = System.getProperty("params");
		String inFasta = System.getProperty("in");
		String inCustomMatrix = System.getProperty("inCustom");
		String distanceName = System.getProperty("distanceName");
		String doPdf = System.getProperty("pdf");
		String zoom = System.getProperty("zoom");
		String UIScaling = System.getProperty("UIScaling");
		String width = System.getProperty("width");
		String height = System.getProperty("height");
		String reference = System.getProperty("reference");
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
				VIEW_SCALE_USER=max(1,new Float(zoom));
			}
			if (UIScaling != null) {
				this.UIScaling = max(1, new Float(UIScaling));
			}
			if (width != null) {
				VIEW_WIDTH = max(1,new Integer(width));
			}
			if (height != null) {
				VIEW_HEIGHT = max(1,new Integer(height));
			}
			if (inCustomMatrix != null) {
				if (distanceName == null) {
					batchError("-DdistanceName must be specified if -DinCustom used.");
				}
				this.distanceName = distanceName;
			}
			//Ok, do it.
			runScript(new Integer(numRuns),new Integer(numIters),params,outputDir,filePref,inFasta,inCustomMatrix,doPdf!=null,reference);
		} else {
			//If we get here, we didn't input any command line arguments.
			System.out.println("Dgraph can be automated: try adding -Dhelp to the arguments when running the program.");
		}
	}
	private void runScript(int numTimes, int numIters, String params1, String outputDirectory, String filePrefix, String inFasta, String customMatrix, boolean pdfOutput, String reference){
		headless = true; //disable rendering while running simulation in script

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
		RUN_MODE = 0;
		if (customMatrix!=null){
			File inputCustomMatrix = new File(customMatrix);
			if (!inputCustomMatrix.exists()){
				batchError("-DinCustom referred to bad/no file\n"+inputCustomMatrix);
			}
			customLoc = customMatrix;
			INPUT_MODE = 2;
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
		needsSetupSimulation = false;
		setupSimulation();

		if (pdfOutput){
			renderUrl = outputDir.toString()+File.separator+filePrefix;
			renderThisFrame = true;
			scriptRenderFrames = 2;
			headless = false;
			return;
		}

		System.out.printf("Setup completed. Running simulations on %d nodes.",nodes.length);
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
			System.out.print(".");
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
			float[][] referenceCoords = new float[numNodes][2];
			if (reference == null || reference.equals("")) {
				getCoordsFile(referenceCoords, filesWrote[bestRunIndex]+".fasta.coords");
			} else {
				getCoordsFile(referenceCoords, reference);
			}
			float[][] otherRun = new float[numNodes][2];
			float[][] out = new float[numNodes][2];
			ana.println("Start of rotations:");
			for(int run = 0; run < numTimes; run++){
				if (reference == null || reference.equals("")) {
					//Rotate all other runs to match run[0]
					if (run==0) continue;
				}
				getCoordsFile(otherRun, filesWrote[run]+".fasta.coords");
				//So, modify filesWrote[run] so that it's better.
				ana.printf("%10d",run);
				PointCloudUtil.transformToMatch2D(referenceCoords, otherRun, out,  ana);
				//Ok. out holds the RMSD minimized transformed otherRun, write it out:
				for(int k = 0; k < nodes.length; k++){
					nodes[k].pos = out[k];
				}
				saveNodeStates(filesWrote[run]+".fasta");
			}
			ana.close();
			
			//Make an alias for best after rotation and render it:
			getCoordsFile(out, filesWrote[bestRunIndex]+".fasta.coords");
			for(int k = 0; k < nodes.length; k++){
				nodes[k].pos = out[k];
			}
			saveNodeStates(outputDir.toString()+File.separator+filePrefix+"-best.fasta");
			renderUrl = outputDir.toString()+File.separator+filePrefix+"-best";
			renderThisFrame = true;
			scriptRenderFrames = 2;
			headless = false;
			return;
		} 
		catch (Throwable e){
			e.printStackTrace();
			batchError("Error in analysis:"+e.toString());
			System.exit(1);
		}
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
		System.out.printf(formatString,"-Dparams","6 numeric simulation parameters seperated by commas:\r\n\t Comparison cutoff\r\n\t Show line cutoff\r\n\t Island removal cutoff \r\n\t Timestep \r\n\t Fluid friction\r\n\t Mass");
		System.out.println();
		System.out.printf(formatString,"-Din","The fasta file to read");
		System.out.println();
		System.out.printf(formatString,"-DinCustom","An NxN matrix (in a text file) of distances to use");
		System.out.println();
		System.out.printf(formatString,"-DdistanceName","Name of the distance metric (required if inCustom set)");
		System.out.println();
		System.out.printf(formatString,"-Dpdf","Flag, means to open the fasta, render a pdf in the same directory, and System.exit(0).");
		System.out.println();
		System.out.printf(formatString,"-Dzoom","Zooming factor. Higher numbers means more zoomed in. Must be >=1.");
		System.out.println();
		System.out.printf(formatString,"-DUIScaling","UI scaling. Bigger numbers mean larger text, etc. Must be >= 1.");
		System.out.println();
		System.out.printf(formatString,"-Dwidth, -Dheight", "width and height of figure in pixels (defaults to 800x800)");
		System.out.println();
		System.out.printf(formatString,"-Dreference", "rotate and flip the final state to best match reference coords file (defaults to match against first run)");
		System.out.println();
	}

	public static final int SHOW_LABELS_SEQUENCE = 0, 
			SHOW_LABELS_NAME = SHOW_LABELS_SEQUENCE + 1, 
			SHOW_LABELS_NUMBER = SHOW_LABELS_NAME + 1, 
			SHOW_LABELS_NONE = SHOW_LABELS_NUMBER + 1;

	private String distanceName = "";
	public String pickUserFile(String title) {
		FileDialog fd = new FileDialog(frame, title, FileDialog.LOAD);
		showFDRememberDirectory(fd);
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
		int windowSize = 0;
		if (INPUT_MODE == -1) {
			String [] optionsStr = { 
					"Fasta sequences, generate PD scores",
					"Fasta sequences + pairwise distance list",
					"Fasta sequences + custom distance matrix",
					"Just labels (one per line) + pairwise distance matrix",
					"Just labels (one per line) + custom distance matrix",
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
			
			SpinnerModel model = new SpinnerNumberModel(22, 7, 100, 1);     
			JSpinner spinner = new JSpinner(model);
			JLabel label = new JLabel("PD score averaging window size (longer than sequence OK):");
			label.setLabelFor(spinner);
			pane.add(label);
			pane.add(spinner);

			JCheckBox alignedSequencePDCB = new JCheckBox("PD score for aligned sequences (overrides windowing)");
			pane.add(alignedSequencePDCB);
			
			options[0].setSelected(true);
			useColoringFileCB.setSelected(false);
			alignedSequencePDCB.setSelected(false);

			boolean shouldRun = JOptionPane.showOptionDialog(
					frame,
					pane,
					"Select input data type",
					JOptionPane.OK_CANCEL_OPTION,
					JOptionPane.PLAIN_MESSAGE,null,null,null) == JOptionPane.OK_OPTION;
			if (!shouldRun) {
				System.exit(00);
			}

			INPUT_MODE = 0;
			for(int k = 0; k < options.length; k++){
				if (options[k].isSelected()){
					INPUT_MODE = k;
					break;
				}
			}
			useColoringFile = useColoringFileCB.isSelected();
			windowSize = alignedSequencePDCB.isSelected() ? 0 : (Integer) spinner.getValue();
		}

		if (INPUT_MODE == 0 || INPUT_MODE == 1 || INPUT_MODE == 2) {
			//Read in a fasta
			if (toRead == null) {
				toRead = pickUserFile("Open a Fasta File");
			}
			inputFileName = new File(toRead).getName();
			String[][] fastaInfo = loadFasta(toRead);
			initNodesFromFasta(fastaInfo);
			if (INPUT_MODE == 0) {
				//Generate distances from fasta
				initDistancesWithPDScore(fastaInfo, windowSize);
			}
		}
		if (INPUT_MODE == 3 || INPUT_MODE == 4) {
			//Custom labels not yet implemented
			throw new RuntimeException("Not yet implemented.");
		}
		if (INPUT_MODE == 1 || INPUT_MODE == 3) {
			//Load pairwise distance list. Fill in the symmetric differences in nodes.
			if (customLoc == null) {
				customLoc = pickUserFile("Select " + distanceName + " matrix file");
				distanceName = JOptionPane.showInputDialog(frame, "Enter the name of your custom score function","Jalview");
			}
			loadPairwiseDistances(customLoc);
		}
		if (INPUT_MODE == 2 || INPUT_MODE == 4) {
			if (customLoc == null) {
				customLoc = pickUserFile("Select " + distanceName + " matrix file");
				distanceName = JOptionPane.showInputDialog(frame, "Enter the name of your custom distance metric","ClustalW");
			}
			loadCustomDistanceMatrix(customLoc);
		}
		if (!headless) {
			//Make a directory in the same directory as toRead to hold interactive screenshots
			renderDir = toRead + "-" + distanceName + "-DGraph-figures";
			new File(renderDir).mkdir();
		}
		if (useColoringFile) {
			setupShaded();
		}
		removeIsolatedNodes();
		//TODO should be a single button to print PDF + PNG.

		//if (!iS) g.beginDraw();
		requestFocus();
	}

	private void loadCustomDistanceMatrix(String customLoc) {
		try {
			Scanner in = new Scanner(openStream(customLoc));
			for(int i = 0; i < nodes.length; i++) {
				nodes[i].distances = new float[nodes.length];
				String[] tokens = in.nextLine().trim().split("\\s+");
				if (tokens.length != nodes.length) {
					throw new RuntimeException("Unexpected line length: " + tokens.length + ". Wanted " + nodes.length + " distances.");
				}
				for(int j = 0; j < nodes.length; j++) {
					nodes[i].distances[j] = new Float(tokens[j]);
				}
			}
			/*
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
			*/
		} catch (Throwable e){
			e.printStackTrace();
			JOptionPane.showMessageDialog(frame, "Error reading distance matrix: " + e.toString(), "Error!", JOptionPane.ERROR_MESSAGE);
			System.exit(1);
		}
	}
	private void loadPairwiseDistances(String customLoc) {
		try {
			//Make a lookup table of node identifiers
			HashMap<String, Integer> nodeLookup = new HashMap<String, Integer>();
			for(int i = 0; i < nodes.length; i++) {
				String fastaLabel = nodes[i].FastaLabel;
				nodes[i].distances = new float[nodes.length];
				nodeLookup.put(fastaLabel, i);
			}
			boolean[][] hasDistance = new boolean[nodes.length][];
			for(int i = 0; i < nodes.length; i++) {
				hasDistance[i] = new boolean[nodes.length];
			}
			Scanner in = new Scanner(openStream(customLoc));
			while(in.hasNextLine()) {
				String a = in.next();
				String b = in.next();
				float score = (float) in.nextDouble(); //not the same as in.nextFloat().
				in.nextLine();
				int idxA = nodeLookup.get(a);
				int idxB = nodeLookup.get(b);
				if (idxA == idxB) {
					continue; //Just ignore pairwise between identical sequences.
				}
				nodes[idxA].distances[idxB] = score;
				nodes[idxB].distances[idxA] = score;
				hasDistance[idxA][idxB] = true;
				hasDistance[idxB][idxA] = true;
			}
			//Check that we got all distances
			for(int i = 0; i < nodes.length; i++) {
				for(int j = 0; j < nodes.length; j++) {
					if (i == j) continue;
					if (!hasDistance[i][j]) {
						throw new RuntimeException("Missing distance between " + nodes[i].FastaLabel + ", " + nodes[i].FastaLabel);
					}
				}
			}
		} catch (Throwable e){
			e.printStackTrace();
			JOptionPane.showMessageDialog(frame, "Error reading pairwise distance list: " + e.toString(), "Error!", JOptionPane.ERROR_MESSAGE);
			System.exit(1);
		}
	}
	//0 window size for computing PDscore with aligned sequences
	private void initDistancesWithPDScore(final String[][] fastaInfo, final int windowSize) {
		distanceName = "PD" + (windowSize == 0 ? "" : " " + windowSize + "-window");
		
		//Compute in parallel
	    int threads = Runtime.getRuntime().availableProcessors();
	    ExecutorService service = Executors.newFixedThreadPool(threads);

		final long[] lastDraw = new long[] {0};
		final int totalDistances = nodes.length * (nodes.length - 1) / 2;
		final int[] distancesComplete = new int[] {0};
		
	    for (int _i = 0; _i < nodes.length; _i++) {
	    	final Integer i = _i;
	        Runnable callable = new Runnable() {
	            public void run() {
	            	int k = i;

	    			float[] dist = new float[nodes.length];
	    			for (int p = 0; p < nodes.length; p++) {
	    				if ((System.nanoTime()-lastDraw[0])>.5e9f) {
	    					lastDraw[0] = System.nanoTime();
	    					if (headless) {
	    						System.err.println("Calculating PDscores: "+nf((distancesComplete[0]/(float)totalDistances)*100, 0, 2)+"%");
		    				} else {
		    					g.beginDraw();
		    					background(0);
		    					fill(255);
		    					textFont(txt);
		    					textAlign(LEFT, TOP);
		    					textSize(txtFontSizeMedium * UIScaling);
		    					text("Calculating PDscores: "+nf((distancesComplete[0]/(float)totalDistances)*100, 0, 2)+"%", 5, 5);
		    					g.endDraw();
		    					repaint();
	    					}
	    				}
	    				
	    				if (p==k) continue;
	    				if (p < k) {
	    					//Symmetry:
	    					continue;
	    				}
	    				distancesComplete[0]++;
	    				if (windowSize == 0) {
	    					if (fastaInfo[1][p].length() != fastaInfo[1][k].length()) {
	    						System.err.println("WARNING: PDscore computation assumed aligned inputs, but " + fastaInfo[0][p] +", " +fastaInfo[0][k] +" have different lengths.");
	    					
	    					}
	    					dist[p]=(float)PDScore.PDScore(fastaInfo[1][p], fastaInfo[1][k]);
	    				} else {
	    					if (fastaInfo[1][k].contains("-")) {
	    						System.err.println("Warning: Sequence " + fastaInfo[0][k] + " contains dashes, windows PDscore is best used on raw sequence (no dashes)");
	    					}
	    					dist[p]=(float)PDScore.PDScoreWindowed(fastaInfo[1][p], fastaInfo[1][k], windowSize);
	    				}
	    			}
	    			nodes[k].distances = dist;
	            }
	        };
	        service.execute(callable);
	    }

	    service.shutdown();
	    try {
	    	service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	    } catch (InterruptedException e) {
	    	System.exit(1);
	    }
	    //Handle symmetry:
	    for(int k = 0; k < nodes.length; k++) {
	    	for(int p = 0; p < k; p++) {
	    		nodes[k].distances[p] = nodes[p].distances[k];
	    	}
	    }
	}
	private void removeIsolatedNodes() {
		removed = "";
		float dontShowCutoff = REMOVAL_ISLAND_FINDER;
		for (int k = 0; k < nodes.length; k++) {
			boolean hasLessThanCutoff = false;
			float[] dist = nodes[k].distances;
			for (int p = 0; p < nodes.length; p++) {
				if (p==k) continue;
				if (dist[p] < dontShowCutoff /**&& dist[p]!=0**/)
					hasLessThanCutoff = true;
			}
			if (!hasLessThanCutoff) {
				removed+=" ,"+(k+1);
				nodes[k].display = false;
			} else {
				nodes[k].display = true;
			}
		}
	}

	private void initNodesFromFasta(String[][] fastaInfo) {
		nodes = new Node[fastaInfo[0].length];
		for (int k = 0; k < nodes.length; k++) {
			nodes[k] = new Node();
		}

		//if (!iS) noLoop();
		for (int k = 0; k < nodes.length; k++) {
			nodes[k].FastaLabel = fastaInfo[0][k];
			nodes[k].FastaSequence = fastaInfo[1][k];
			nodes[k].myId = k;
			final String hashCoordText = "#COORDS:";
			int hashCoords = nodes[k].FastaLabel.indexOf(hashCoordText);
			if (hashCoords==-1) { //No coords specified, pick a random location
				scatterNode(nodes[k]);
			} 
			else { //Coords specified, strip from the label
				String[] pos = nodes[k].FastaLabel.substring(hashCoords+hashCoordText.length()).split(","); //Split on comma for x,y
				nodes[k].pos[0] = new Float(pos[0]);
				nodes[k].pos[1] = new Float(pos[1]); 
				nodes[k].FastaLabel = nodes[k].FastaLabel.substring(0, hashCoords).trim(); //The "Real" header, remove trailing leading whitespace.
				//System.out.printf("Initial position: %s (%.3f,%.3f)",nodes[k].FastaLabel,nodes[k].pos[0],nodes[k].pos[1]);
				//System.out.println();
			}
		}
	}

	private String[][] _loadFasta(String[] fastaLines) {
		String[][] fastaInfo = null;
		try {
			fastaInfo = unFastaSuper(fastaLines);
		} catch (Throwable e) {
			e.printStackTrace();
			JOptionPane.showMessageDialog(frame, "Error reading fasta file: " + e.toString(), "Error!", JOptionPane.ERROR_MESSAGE);
			System.exit(1);
		}
		return fastaInfo;
	}
	
	private String[][] loadFasta(String toRead) {
		String[] lines = loadStrings(toRead);
		return _loadFasta(lines);
	}

	public void scatterNode(Node got) {
		float theta = random(TWO_PI);
		float mag = random(0.01f, 1);
		got.pos[0] = mag*cos(theta);
		got.pos[1] = mag*sin(theta);
		got.vel[0] = 0;
		got.vel[1] = 0;
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
		float amtBlue = pow(constrain(-2*var+1+.4f, 0, 1), 2);
		float amtCream = sqrt(constrain(-1.5f*(1-var)+1, 0, 1));
		int c = lerpColor(lerpColor(color(100, 100, 255) /*light blue*/, color(100, 170, 100) /*faded green*/, (1-amtBlue)), color(255, 255, 152) /*cream*/, amtCream);
		stroke(lerpColor(c, color(255, 255, 255) /*white*/, var*.5f));
	}
	private Node[] nodes;
	private float TOTAL_VARIANCE;
	private float RANDOM_VARIANCE = -1;
	private final float ZOOM = 50; //ZOOM effects the scaling of the position integration; it used to be determined by screen coords ...
	private float VIEW_SCALE = -1; //Updated from VIEW_SCALE_USER
	private float VIEW_SCALE_USER = 25;
	private int showLabels = SHOW_LABELS_NUMBER; 
	private boolean showLines = true;
	private int NodeDotSize = 3; //rectangle of this radius at each node
	private String underMouse;
	int simulationTicks = 0;

	private boolean UPDATES_ON_DRAW = true; //Normally true, but set to false on the first frame in the GUI.
	private boolean wantsFrameBreak = false;
	private boolean simStarted = true;
	public void drawSimulation(boolean controlInfo) {
		background(255);
		UPDATES_ON_DRAW = simStarted;
		underMouse = "";
		textFont(txt);
		textSize(txtFontSizeMedium * UIScaling);
		//Start out with some sane defaults:
		noStroke();
		fill(0);
		stroke(0);
		float padding = 5;
		StringWriter sw;
		PrintWriter pw;
		//New textbox
		sw = new StringWriter(); pw = new PrintWriter(sw);
		pw.println("Input file: " + inputFileName);
		pw.println("\"Islands\" Removed by "+REMOVAL_ISLAND_FINDER+" distance threshold: "+removed.substring(min(removed.length(), 2)));
		//2, controlInfo?20:5, width, 80);
		//pw.println("Noncomparable Seperation Factor:"+CURRENT_CUTOFF); //, 2, 4);
		//pw.println("Displaying:"+SHOW_LINES); //, 2, 100);
		textAlign(LEFT, TOP);
		fill(0);
		text(sw.toString(), padding, padding);
		
		if (controlInfo) {
			if (!UPDATES_ON_DRAW) {
				fill(255, 0, 0);
				textAlign(RIGHT, TOP);
				text("Simulation stopped. Press 'y' to start again.", width - padding, padding);
			}
			if (SHOW_HELP_TEXT) {
				sw = new StringWriter(); pw = new PrintWriter(sw);
				pw.println("'h' to hide/show this help");
				pw.println("'s' save figure (.PDF & .PNG)");
				pw.println("'m' save figure (.COORDS)");
				pw.println("'l' to switch on/off lines");
				pw.println("'c' to turn off labels");
				pw.println("'n' to show sequence numbering");
				pw.println("'a' to show sequences");
				pw.println("'b' to show sequence headers");
				pw.println("'9','0' zooms out, in");
				pw.println("'7','8' makes text bigger, smaller");
				pw.println("'j' sets the physics options");
				pw.println("'d' delete nodes under mouse");
				pw.println("'q', 'w' make dots bigger, smaller");
				pw.println("'e' to scatters the points.");
				pw.println("'y' to pause/resume.");

				fill(0);
				textAlign(RIGHT, BOTTOM);
				text(sw.toString(), width - padding, height - padding);
			}
		}
		if (true) {
			fill(0);
			textAlign(LEFT, BOTTOM);
			text(distanceName+" distance", 20, height * 0.25f - (textAscent() + textDescent())/2 - padding);
			beginShape(QUADS);
			float minY = height * 0.25f;
			float maxY = height * 0.75f;
			for (float PD = 0; PD <= SHOW_LINES; PD+=SHOW_LINES/20) {
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
			for (float PD = 0; PD <= SHOW_LINES; PD+=(SHOW_LINES/20) * 5) {
				float mag = PD/SHOW_LINES;
				float y = lerp(minY, maxY, mag);
				textAlign(LEFT, CENTER);
				text(PD, 30, y);
			}
		}
		//Draw simulation:
		VIEW_SCALE = VIEW_SCALE_USER * VIEW_SCALE_USER;
		//keyPressed2();
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
		translate(-width/2, -height/2);
		
		if (controlInfo) {
			fill(0);
			textSize(txtFontSizeMedium * UIScaling);
			textAlign(LEFT, BOTTOM);
			text("Simulation iterations:"+simulationTicks+"\nSolution score (lower is better):"+TOTAL_VARIANCE+"\nTimes smaller than random:"+RANDOM_VARIANCE/TOTAL_VARIANCE+"\n"+
					"Zoom: " + VIEW_SCALE_USER + ", UIScaling: " + UIScaling + "\n" +
					"Nodes under mouse:"+underMouse.substring(min(underMouse.length(), 2)), padding, height - padding);
			stroke(0);
		}
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
		loop();
		if (keyPressed && Character.toLowerCase(key)=='y') {
			acted = System.nanoTime();
			simStarted = !simStarted;
			return;
		}
		if (keyPressed && Character.toLowerCase(key)=='s') {
			acted = System.nanoTime();
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
		/*
		if (keyPressed && Character.toLowerCase(key)=='r') {
			acted = System.nanoTime();
			setupSimulation();
			return;
		}
		*/
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


		if ((keyPressed && Character.toLowerCase(key)=='n')) {
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
			SHOW_HELP_TEXT = !SHOW_HELP_TEXT;
			return;
		}

		if ((keyPressed && Character.toLowerCase(key)=='l')) {
			acted = System.nanoTime();
			showLines = !showLines;
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
			VIEW_SCALE_USER=max(VIEW_SCALE_USER*0.95f, 1);
			return;
		}


		if ((keyPressed && Character.toLowerCase(key)=='0')) {
			acted = System.nanoTime();
			VIEW_SCALE_USER=VIEW_SCALE_USER*1.05f;
			return;
		}

		if ((keyPressed && Character.toLowerCase(key)=='7')) {
			acted = System.nanoTime();
			UIScaling=max(UIScaling - 0.05f,0.25f);
			return;
		}

		if ((keyPressed && Character.toLowerCase(key)=='8')) {
			acted = System.nanoTime();
			UIScaling=min(UIScaling + 0.05f,4f);
			return;
		}

		if ((keyPressed && Character.toLowerCase(key)=='j')) {
			acted = System.nanoTime();
			new Thread() {
				public void run() {
					newPhysicsVariables();
				}
			}.start();
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
	float UIScaling = 1.0f;

	/** DISTANCE VARIABLES **/

	private float CURRENT_CUTOFF = 14f;
	private float SHOW_LINES = 14f;
	private float REMOVAL_ISLAND_FINDER = 14f;

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
			float newRemovalIslandFinder = PApplet.parseFloat(areas[2].getText());
			boolean needsNewIsolatedNodes = newRemovalIslandFinder != REMOVAL_ISLAND_FINDER;
			REMOVAL_ISLAND_FINDER = newRemovalIslandFinder;
			if (needsNewIsolatedNodes) {
				removeIsolatedNodes();
			}
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
		boolean display = true;
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
			translate(offx, offy);


			cx = screenX(0, 0);
			cy = screenY(0, 0);
			if (pow(mouseX - cx, 2)+pow(mouseY-cy, 2) < 3000) {
				translate(-offx, -offy);
				underMouse+=", "+(myId+1);
				return;
			}

			//Label will be placed distance r away, at some angle, in a circle of radius circleSize
			float circleSize = 10*UIScaling;
			float r = max(25*UIScaling, NodeDotSize + 5);
			//Try a random sweep around the circle 
			int anglesToTest = 20;
			//Default to an angle that's away from the center of the picture (which is at 0,0)
			float angle0 = atan2(offy, offx);
			for(int i = 0; i <= anglesToTest; i++) {
				boolean canDraw = true;
				float angle;
				if (i == anglesToTest) {
					//Last iteration. Use the default angle, ignore overlap.
					//angle = angle0;
					//Last iteration, use a seeded random angle:
					//angle = (labelString.hashCode() % TWO_PI);
					angle = angle0 + labelString.hashCode() % (TWO_PI / anglesToTest);
				} else {
					//Flip back and forth around angle0 at intervals:
					angle = angle0 + (i % 2 == 0 ? - 1 : 1) * ((i + 1)&(~1)) * TWO_PI / anglesToTest;
					//Angles close to horizontal are ugly, and lead to unreadable labels when showing sequence
					//as these overflow the "circleR" space at the label. So avoid that:
					float horizThreshold = TWO_PI/20;
					if (abs(angle + horizThreshold) % PI < horizThreshold * 2) {
						canDraw = false;
					}
					//Check drawability slightly to the left and right of the proposed label position
					if (!labelPlanner.canDrawAt(cx + r*cos(angle) - circleSize, cy + r*sin(angle))) {
						canDraw = false;
					}
					if (!labelPlanner.canDrawAt(cx + r*cos(angle) + circleSize, cy + r*sin(angle))) {
						canDraw = false;
					}
					//Check drawability on the segment connecting label to point. Can't check at the center of the ray since we want to allow
					//labels for points radiating from the exact same spot.
					for (float rad = max(10, NodeDotSize + 5); rad < r + circleSize && canDraw; rad+= 5) {
						float ox = rad*cos(angle);
						float oy = rad*sin(angle);
						if (!labelPlanner.canDrawAt(cx + ox, cy + oy)) {
							canDraw = false;
						}
					}
				}
				if (canDraw) {
					float ox = r*cos(angle);
					float oy = r*sin(angle);

					stroke(0, 0, 0, 100);
					line(0, 0, ox*((r-circleSize)/r), oy*((r-circleSize)/r));
					fill(1, 1, 1, 40);
					if (showLabels == SHOW_LABELS_NUMBER) {
						ellipse((int)ox, (int)oy, circleSize*2, circleSize*2);
					}
					fill(0);
					translate(-offx, -offy);
					textFont(txt);
					textSize(txtFontSizeSmall * UIScaling);
					textAlign(LEFT, BOTTOM);
					text(labelString, (int)(offx+ox-textWidth(labelString)/2), (int)(offy+oy + textAscent()/2 + textDescent()));
					labelPlanner.mark(cx + ox, cy + oy, circleSize);
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
		public void drawLineTo(int oth, float goodness) {
			//called by update, but if we're in script, don't do it!
			if (headless) return;
			//Done.
			boolean trans = transparentCountdown>0;
			if (!trans) {
				if (showLines && distances[oth]<SHOW_LINES) {
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
					continue;
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




	public void runUtility(int whichUtility){
		headless = true;

		FileDialog fd = new FileDialog(frame, "Please select your input file", FileDialog.LOAD);
		showFDRememberDirectory(fd);
		File got = null;
		String toRead = "";
		if ( fd.getDirectory() != null && fd.getFile() != null) 
		{
			try {
				got = new File( fd.getDirectory() + File.separator + fd.getFile());
				if (!got.exists()){
					JOptionPane.showMessageDialog(frame, "You entered a nonexistant file:\n"+got);
					runUtility(whichUtility);
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
		
		String[] utilitySuffixes = new String[] {
				".fasta",
				"-sequences.txt",
				"-unique.txt",
				"-headers.txt",
				"-distsfromclustalw.txt",
				"-distsfrompdunaligned.txt",
				"-distsfrompdaligned.txt",
				"-distsfromdnasimilarity.txt",
				"-sequencesislandsremoved.txt",
				"-distsfromjalviewpairwisealn.txt",
				"-reordered.fasta",
				"-seqnumbers.txt",
		};

		// String gotString = JOptionPane.ShowInputDialog("Please input a meaningful description of this graph.");
		String[] lines = loadStrings(toRead);
		File outtt = new File(got.toString()+utilitySuffixes[whichUtility]);
		try {
			outtt.createNewFile();
			System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(outtt))));
		} 
		catch (Throwable e){
			e.printStackTrace();
			die("");
		};
		String errorMessage = "";
		try {
			if (whichUtility==0){
				Fastize(lines); 
				return;
			}
			if (whichUtility==1){
				unFasta(lines,false,false); 
				return;
			}
			if (whichUtility==2){
				removeDuplicates(lines); 
				return;
			}
			if (whichUtility==3){
				unFasta(lines,true,false); 
				return;
			}
			if (whichUtility==4){
				ClustalwUtil.clustalwAnalyze(lines); 
				return;
			}
			if (whichUtility==5){
				pdScore2Mat(lines);
				return;
			}

			if (whichUtility==6){
				pdScoreMat(lines);
				return;
			}

			if (whichUtility==7){
				seqSimScore(lines);
				return;
			}

			if (whichUtility==8){
				pdScoreMultisearch(lines);
				return;
			}

			if (whichUtility==9){
				JalviewUtil.jalviewPairwiseAlignmentScores(lines); 
				return;
			}

			if (whichUtility==10){
				fd = new FileDialog(frame, "Select a fasta as the ordering reference", FileDialog.LOAD);
				showFDRememberDirectory(fd);
				String[] ordering_reference;
				if ( fd.getDirectory() != null && fd.getFile() != null) 
				{
					got = new File( fd.getDirectory() + File.separator + fd.getFile());
					if (!got.exists()){
						JOptionPane.showMessageDialog(frame, "You entered a nonexistant file:\n"+got);
						System.exit(1);
					}
					ordering_reference = loadStrings(got.getAbsoluteFile().toString());
				}  
				else {
					//We canceled the FD - exit normally.
					System.exit(0);
					//Blech.
					return;
				}
				
				FastaUtil.reorderFasta(lines, ordering_reference); 
				return;
			}
			if (whichUtility==11){
				unFasta(lines,true,true); 
				return;
			}
		} catch (Throwable e) {
			errorMessage = "Errors occurred: " + e.toString();
			throw new RuntimeException(e);
		} finally {
			System.out.flush();
			JOptionPane.showMessageDialog(frame,"O.K. " + errorMessage + " the output has been stored to "+outtt);
		}
	}

	private String fdLastDirectory = System.getProperty("user.home");
	private void showFDRememberDirectory(FileDialog fd) {
		fd.setDirectory(fdLastDirectory);
		fd.show();
		String chosenDirectory = fd.getDirectory();
		if (chosenDirectory != null) {
			fdLastDirectory = chosenDirectory;
		}
	}

	public void pdScore2Mat(String[] lines){
		int wSize = -1;
		String input =  JOptionPane.showInputDialog(frame,"Input the Window Size(integer):","22");
		try {
			wSize = new Integer(input);
		} 
		catch (Throwable e){
			JOptionPane.showMessageDialog(frame,("Non integer input."),"You entered a non-integer window size.",JOptionPane.ERROR_MESSAGE);
			return;
		}
		String[][] fastaInfo = _loadFasta(lines);
		initNodesFromFasta(fastaInfo);
		initDistancesWithPDScore(fastaInfo, wSize);
		for(int k = 0; k < nodes.length; k++){
			for(int p = 0; p < nodes.length; p++){
				System.out.print(nodes[k].distances[p]+(p < nodes.length - 1 ? "\t" : ""));
			}
			System.out.println();
		}
	}

	public void pdScoreMat(String[] lines){
		String[][] fastaInfo = _loadFasta(lines);
		initNodesFromFasta(fastaInfo);
		initDistancesWithPDScore(fastaInfo, 0);
		for(int k = 0; k < nodes.length; k++){
			for(int p = 0; p < nodes.length; p++){
				System.out.print(nodes[k].distances[p]+(p < nodes.length - 1 ? "\t" : ""));
			}
			System.out.println();
		}
	}

	public void pdScoreMultisearch(String[] lines){
		FileDialog fd = new FileDialog(frame, "Choose the short list of peptides (to match against)", FileDialog.LOAD);
		showFDRememberDirectory(fd);
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
	public void unFasta(String[] lines, boolean printHeaders, boolean printNumbers){
		ArrayList unFasta = new ArrayList();
		ArrayList headers = new ArrayList();
		unFasta0(unFasta,headers,lines);
		for(int k = 0; k < unFasta.size(); k++){
			if (printNumbers) {
				System.out.printf("%-8d ", k+1);
			}
			if (printHeaders){
				System.out.println(headers.get(k));
			} 
			else {
				System.out.println(unFasta.get(k));
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
				showFDRememberDirectory(fd);
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

			try {
				PrintWriter out = new PrintWriter(new FileWriter(got));
				for(int k = 0; k < nodes.length; k++){
					if (!nodes[k].display) continue;
					out.println(">"+nodes[k].FastaLabel+" #COORDS:"+nodes[k].pos[0]+","+nodes[k].pos[1]);
					String seq = nodes[k].FastaSequence;
					for(int c = 0; c < seq.length(); c+=80){
						out.println(seq.substring(c,min(seq.length(),c+80)));
					}
				}
				out.close();

				out = new PrintWriter(new FileWriter(new File(toRead+".coords")));
				for(int k = 0; k < nodes.length; k++){
					if (!nodes[k].display) continue;
					out.printf("%d %f %f %s",k+1,nodes[k].pos[0],nodes[k].pos[1],nodes[k].FastaLabel);
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
			if (!headless){
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
				System.out.printf("%f ",seqSimDist((String)unFasta.get(k),(String)unFasta.get(p)));   
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
			//may not equal k since there can be islands that we don't export
			int whichNode = new Integer(coord[0]) - 1;
			coords[whichNode][0] = new Float(coord[1]);
			coords[whichNode][1] = new Float(coord[2]);
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
