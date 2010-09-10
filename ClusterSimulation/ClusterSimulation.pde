import processing.pdf.*;

String removed = "";


import java.awt.*;

boolean PRESETUP = false;
int WHAT_DO_YOU_WANT = -1;
boolean singleFramePdf = false;

public ClusterSimulation(){  
  doCommandLine();
  if (singleFramePdf){
    PRESETUP = true;
    WHICH_SCREEN = 2;
  }
  if (!PRESETUP){
    PRESETUP = true;
    //Show a dialog with options of what to do:
    String [] optionsStr = { 
      "The Simulation Major < YOU PROBABLY WANT THIS!", 
      "Utility: Sequence list -> Fasta Format", 
      "Utility: Fasta Format -> Sequence List", 
      "Utility: Remove Duplicates from Sequence List", 
      "Utility: Fasta Format -> List of 'Fasta Headers'",
      "Utility: Clustalw output -> 'distance' matrix (experimental!)",
      "Utility: Sequence List -> Neo-PD Score Matrix (Long Proteins)",
      "Utility: Sequence List -> PD Score Matrix (Short Peptide)",
      "Utility: DNA Fasta -> Similarity Score Matrix",
      "Multisearch: Choose from a list peptides similar to any in another list",
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

    JOptionPane.showMessageDialog(
    holder,
    pane,
    "Welcome! What would you like to run?",
    JOptionPane.INFORMATION_MESSAGE,null);

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
}
void setup(){
  
  //Test here if we have write priviledges
  String test = sketchPath+File.separator+".test";
  try {
    FileOutputStream fos = new FileOutputStream(test);
    fos.write(2);
    fos.close();
  } 
  catch (Throwable e){
    JOptionPane.showMessageDialog(frame,"Please copy this program to a directory that has file-write access.\nNo installation is necessary except for this step.\n"+sketchPath+" does not have write access.","Write-Protected Directory",JOptionPane.ERROR_MESSAGE);
    System.exit(1);
  }
  //The following called twice, if that makes any sense.
  if (WHAT_DO_YOU_WANT==-1){
    if (!frame.isResizable()){
      frame.setResizable(true);
    }
    size(800,800,P2D);
    txt = createFont("Dialog",30,true);
  } 
  else {
    TestFasta(WHAT_DO_YOU_WANT);
    exit();
  }
}

int WHICH_SCREEN = 0;
void draw(){
  frameRate(60);
  if (textPaneClone==null){
    textPaneClone = createGraphics(2048,2048,P2D);
    textPaneClone.beginDraw();
  }
  if (WHICH_SCREEN == 1){
    frame.setTitle("ClusterSimulation: dGraph copyright Ben Braun, 2009");
    UpdateSimulation(true);
    if (keyPressed && keyEvent.getKeyChar()=='p'){
      WHICH_SCREEN++;
    }
  } 
  else if (WHICH_SCREEN==2){
    DrawShaded();
    if (singleFramePdf){
      System.exit(0);
    }
  }
  else {
    if (frameCount <= 3){
      background(254,200,224);
    }
    if (frameCount==3){
      reinit();
      WHICH_SCREEN = 1;
    }
  }
}












