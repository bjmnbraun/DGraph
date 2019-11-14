import java.awt.event.*;
import java.util.*;

public static final int SHOW_LABELS_SEQUENCE = 0, 
SHOW_LABELS_NAME = SHOW_LABELS_SEQUENCE + 1, 
SHOW_LABELS_NUMBER = SHOW_LABELS_NAME + 1, 
SHOW_LABELS_NONE = SHOW_LABELS_NUMBER + 1;

private String toRead = null;
private String scoreName = "";
String pickUserFile(String title) {
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
      return got.getAbsoluteFile().toURI().toURL().toString();
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
void reinit() {
  if (!iS) g.endDraw();

  simulationTicks = 0;
  /*
  if (toRead!=null && !(new File(toRead).exists())){
   toRead = null;
   }
   */
  if (toRead==null) {
    toRead = pickUserFile("Open a Fasta File");
    try {
      Thread.sleep(100);
    } 
    catch (Throwable e) {
    };
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

  if (!iS) g.beginDraw();

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
  if (!iS) noLoop();
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
      float PIX_TO_PD_FACTOR = (ZOOM/640.);
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
    if (!iS) if ((System.nanoTime()-lastDraw)>.5e9) {
      lastDraw = System.nanoTime();
      background(0);
      fill(255);
      textFont(txt);
      textAlign(LEFT, TOP);
      textSize(30);
      text("Calculating PDscores:"+nf((k/(float)nodes.length)*100, 0, 2)+"%", 5, 5);
      g.endDraw();
      repaint();
      g.beginDraw();
    }
  }
  //Phew.
  if (!iS) {
    loop();
    g.endDraw();
    frame.setSize(new Dimension(frame.getWidth()-10, frame.getHeight()-10));
    try {
      Thread.sleep(100);
    } 
    catch (Throwable e) {
    };
    frame.setSize(new Dimension(frame.getWidth()+10, frame.getHeight()+10));
    g.beginDraw();
  }
}
void scatterNode(Node got) {
  float theta = random(TWO_PI);
  float mag = random(640/200f, 640/2f);
  got.pos[0] = mag*cos(theta);
  got.pos[1] = mag*sin(theta);
}
void scatterNodes() {
  //Similar to reinit, except that only the positions of the already existant nodes are moved.
  for (int k = 0; k < nodes.length; k++) {
    Node got = nodes[k];
    if (!got.exists) continue;
    scatterNode(got);
  }
}
void strokeColors(float val) {
  float var = constrain(val/SHOW_LINES, 0, 1);
  float amtRed = pow(constrain(-2*var+1+.4, 0, 1), 2);
  float amtYello = constrain(1-3*abs(var-.3), 0, 1);
  float amtBlue = sqrt(constrain(-1.5*(1-var)+1, 0, 1));
  int c = lerpColor(lerpColor(color(0, 0, 255), color(0, 170, 40), (1-amtRed)), color(255, 255, 152), amtBlue);
  stroke(lerpColor(c, color(255, 255, 255), var*.5));
}
private Node[] nodes;
private float TOTAL_VARIANCE;
private float RANDOM_VARIANCE = -1;
private final float ZOOM = 50; //ZOOM effects the scaling of the position integration; it used to be determined by screen coords ...
private float VIEW_SCALE = -1; //Updated from VIEW_SCALE_USER
private float VIEW_SCALE_USER = 25;
private int showLabels = SHOW_LABELS_NONE; 
private boolean ShowMinLines = false;
private int NodeDotSize = 7; //px
private String underMouse;
int simulationTicks;
PGraphics textPaneClone;

private boolean UPDATES_ON_DRAW = true; //Normally true, but set to false on the first frame in the GUI.
private boolean wantsFrameBreak = false;
private boolean simStarted = false;
void UpdateSimulation(boolean controlInfo) {
  UPDATES_ON_DRAW = simStarted; //GUI uses this so the user can just look at it for a second.
  underMouse = "";
  background(255);
  if (showLabels != SHOW_LABELS_NONE) {
    textPaneClone.background(255);
  }
  textFont(txt);
  textPaneClone.textFont(txt);
  textSize(20);
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
    textSize(14);
    int fromBottomRight = 150;
    text("'h' to switch on/off lines", width-5, height-fromBottomRight-33);
    text("'c' to turn off labels", width-5, height-fromBottomRight-20);
    text("'s' to show sequence numbering", width-5, height-fromBottomRight-7);
    text("'a' to show sequences", width-5, height-fromBottomRight+5);
    text("'b' to show sequence headers", width-5, height-fromBottomRight+18);
    text("'9' zooms out", width-5, height-fromBottomRight+35);
    text("'0' zooms in", width-5, height-fromBottomRight+45);
    textSize(16);
    text("'m' saves the current state", width-5, height-fromBottomRight+55);
    text("'p' proceeds to the analysis step", width-5, height-fromBottomRight+72);
    text("'j' sets the physics options", width-5, height-fromBottomRight+85);
    text("'d' disable nodes", width-5, height-fromBottomRight+100);
    textSize(12);
    text("'w' increases dot-size", width-5, height-fromBottomRight+115);
    text("'q' decreases dot-size", width-5, height-fromBottomRight+125);
    text("'r' or 'e' re-run the simulation. 'e' only scatters the points, 'r' does more.", width-5, height-fromBottomRight+135);
  }
  if (true) {
    if (transparentCountdown>0) {
      textSize(12);
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
      textSize(12);
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
  textPaneClone.pushMatrix();
  textPaneClone.translate(width/2, height/2);
  translate(width/2, height/2);
  //The actual simulation:
  runActualSimulation();
  //End actual simulation.
  flushLines();
  if (RANDOM_VARIANCE<0) {
    RANDOM_VARIANCE = TOTAL_VARIANCE;
  }
  for (int k = 0; k < nodes.length; k++) {
    nodes[k].draw();
  }
  for (int k = 0; k < nodes.length; k++) {
    nodes[k].drawDot();
  }
  wantsFrameBreak = false;
  //popMatrix();
  textPaneClone.popMatrix();
  //image(textPaneClone,0,0,width,height,0,0,width,height);
  transparentCountdown--;
  fill(0);
  textSize(12);
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
}
private long acted = System.nanoTime();
public void keyPressed(java.awt.event.KeyEvent e) {    
  //  super.keyPressed(e);
  //System.out.println("Breaking..."+System.nanoTime());
  wantsFrameBreak = true;
  key = e.getKeyChar();
  keyPressed = true;
  if (System.nanoTime() - acted > .033e9) {
    keyPressed2();
  }
}
public void keyPressed2() {
  if (keyPressed && Character.toLowerCase(key)=='y') {
    acted = System.nanoTime();
    simStarted = !simStarted;
    return;
  }
  if (keyPressed && Character.toLowerCase(key)=='d') {
    acted = System.nanoTime();
    boolean canDo = JOptionPane.showConfirmDialog(frame, "Disable all nodes under your mouse? This cannot be undone.")==JOptionPane.OK_OPTION;
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
    reinit();
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
    VIEW_SCALE_USER=VIEW_SCALE_USER+1;
    return;
  }

  if ((keyPressed && Character.toLowerCase(key)=='0')) {
    acted = System.nanoTime();
    VIEW_SCALE_USER=max(VIEW_SCALE_USER-1, 1);
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
  int compareTo(Object oth) {
    return (int)Math.signum(myz-((LineSort)oth).myz);
  }
}
int lineRenderingUsed = 0;
private LineSort[] lineRendering;
void flushLines() {
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

/** DISTANCE VARIABLES **/

private float CURRENT_CUTOFF = 999;
private float SHOW_LINES = 14;
private float REMOVAL_ISLAND_FINDER = 12;

/** PHYSICS VARIABLES **/

float timeStep = .001;
float fluidFriction = .008;
float mass = .01;

abstract class VarSetter extends JPanel {
  public VarSetter(int numSettings) {
    this.numSettings = numSettings;
    areas = new JTextArea[numSettings];
    labels = new JLabel[numSettings];
  }
  int numSettings;
  JTextArea[] areas;
  JLabel[] labels;
  abstract void dedicate();
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
  void dedicate() {
    timeStep = float(areas[0].getText());
    fluidFriction = constrain(float(areas[1].getText()), 0, 1);
    mass = float(areas[2].getText());
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
  void dedicate() {
    CURRENT_CUTOFF = float(areas[0].getText());
    SHOW_LINES = float(areas[1].getText());
    REMOVAL_ISLAND_FINDER = float(areas[2].getText());
  }
}
void newPhysicsVariables() {
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

/** END PHYRSICS **/

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
  void draw() {
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
    flush();
  }
  void drawDot() {
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
      drawLineTo0(nodeConnection, brightness(shadeColor)/60.);
      nodeConnection = -1;
    }
    rect(x-1, y-1, NodeDotSize, NodeDotSize);
    flush();
  }
  void placeText(String stuff, float offx, float offy) {
    float cx=0, cy=0;
    for (int k = -1; k <= 0; k++) {
      if (k==0) {
        PApplet e = ClusterSimulation.this;
        e.fill(0);
        e.textFont(txt);
        e.textSize(12);
        e.textAlign(CENTER, CENTER);
        e.translate(offx, offy);
      } 
      else {
        PGraphics e = textPaneClone;
        e.fill(0);
        e.textFont(txt);
        e.textSize(12);
        e.textAlign(CENTER, CENTER);
        e.translate(offx, offy);
      }
    }
    for (int k = -1; k <= 0; k++) {
      if (k==0) {
        PApplet e = ClusterSimulation.this;
        cx = screenX(0, 0);
        cy = screenY(0, 0);
        if (pow(mouseX - cx, 2)+pow(mouseY-cy, 2) < 3000) {
          e.translate(-offx, -offy);
          if (k==0) {
            underMouse+=", "+(myId+1);
            return;
          }
        }
      } 
      else {
        PGraphics e = textPaneClone;
        cx = screenX(0, 0);
        cy = screenY(0, 0);
        if (pow(mouseX - cx, 2)+pow(mouseY-cy, 2) < 3000) {
          e.translate(-offx, -offy);
          if (k==0) {
            underMouse+=", "+myId;
            return;
          }
        }
      }
    }

    int white = color(255, 255, 255);
    float r0 = atan2(offy, offx);
    int poss = 10;
    for (float r = r0; r < r0+TWO_PI; ) {
      boolean bad = false;
      for (float rad = 10; rad < 40; rad+= 5) {
        float ox = cx+rad*cos(r);
        float oy = cy+rad*sin(r);
        if (textPaneClone.pixels[constrain((int)ox+((int)oy)*textPaneClone.width, 0, textPaneClone.pixels.length-1)]!=white) {
          bad = true;
        }
      }
      float sidethresh = TWO_PI/poss/2;
      if (abs((r % PI) + PI) < sidethresh || abs((r % PI) - 0) < sidethresh || abs((r%PI) - PI) < sidethresh) {
        bad = true;
      }
      boolean oBad = false && (r+TWO_PI/poss)>TWO_PI;
      if (!bad || oBad) {
        float goRad = !oBad?25:40;
        float ox = goRad*cos(r);
        float oy = goRad*sin(r);
        for (int k = -1; k <= 0; k++) {
          if (k==0) {
            PApplet e = ClusterSimulation.this;
            float circleSize = 20;
            e.stroke(0, 0, 0, 100);
            e.line(0, 0, ox*((goRad-circleSize/2)/goRad), oy*((goRad-circleSize/2)/goRad));
            e.fill(1, 1, 1, 40);
            if (showLabels == SHOW_LABELS_NUMBER) {
              e.ellipse(ox, oy, circleSize, circleSize);
            }
            e.fill(0);
            e.translate(-offx, -offy);
            //float x = pos[0]*cos(pos[1])/(width/(40));
            //float y = pos[0]*sin(pos[1])/(width/(40));
            if (k==0) {
              e.text(stuff/*+" ("+nf(x,0,2)+","+nf(y,0,2)+")"*/, (int)(offx+ox), (int)(offy+oy));
              return;
            }
          }
          else {
            PGraphics e = textPaneClone;
            float circleSize = 20;
            e.stroke(0, 0, 0, 100);
            e.line(0, 0, ox*((goRad-circleSize/2)/goRad), oy*((goRad-circleSize/2)/goRad));
            e.fill(1, 1, 1, 40);
            e.ellipse(ox, oy, circleSize, circleSize);
            e.fill(0);
            e.translate(-offx, -offy);
            //float x = pos[0]*cos(pos[1])/(width/(40));
            //float y = pos[0]*sin(pos[1])/(width/(40));
            //if (k==0){
            e.text(stuff/*+" ("+nf(x,0,2)+","+nf(y,0,2)+")"*/, (int)(offx+ox), (int)(offy+oy));
            //return;
            //}
          }
        }
      }
      r+=TWO_PI/poss;
    }
    /* //Doesn't happen anymore.
     fill(255,255,255,40);
     ellipse(0,0,15,15);
     fill(0);
     text(stuff,0,0);
     */
    for (int k = -1; k <= 0; k++) {
      if (k==0) {
        PApplet e = ClusterSimulation.this;
        e.translate(-offx, -offy);
      }
      else {
        PGraphics e = textPaneClone;
        e.translate(-offx, -offy);
      }
    }
  }
  void setStroke(float goodness, float distance) {
    strokeWeight(2);
    colorMode(HSB, 255);
    stroke(200+goodness/SHOW_LINES*50, 255, 255);
    colorMode(RGB, 255);
  }
  boolean drawsLines() {
    return (transparentCountdown>0 || ShowMinLines) && !wantsFrameBreak;
  }
  void drawLineTo(int oth, float goodness) {
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
  void drawLineTo0(int oth, float goodness) {
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
  void update() {
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
      float proximity = (compLength*(ZOOM/640.)); //pixel to pd
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

      if (proximity < 1e-5) { //divide by zero catch
        dispAdd[0] = -dispar*1e-2;
        dispAdd[1] = -dispar*1e-2;
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
  float[] componentsTo(float[] other) {
    componentsMemSaver[0] = other[0]-pos[0];
    componentsMemSaver[1] = other[1]-pos[1];
    return componentsMemSaver;
  }
}
float[] dispAddMemorySaver = new float[2];
float[] componentsMemSaver = new float[2];
void addVec(float[] a, float[] b) {
  a[0]+=b[0];
  a[1]+=b[1];
}
void addVec(float[] a, float b1, float b2) {
  a[0]+=b1;
  a[1]+=b2;
}
