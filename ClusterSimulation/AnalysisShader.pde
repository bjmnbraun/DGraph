import javax.swing.*;

float[] shadings;
int[] nodeConnections;
long gotHere;
String pdfUrl = null;
void DrawShaded(){
  boolean renderThisFrame = false;
  //Precalc
  String outputPdfDir = null;
  if (shadings==null){
    renderThisFrame = true;
    gotHere = System.nanoTime();
    shadings = new float[nodes.length]; //defaults to 0
    nodeConnections = new int[nodes.length];

    String toRead = null;
    String[] commands = null;
    if (singleFramePdf){
      commands = new String[0];
    }
    if (commands==null){
      FileDialog fd = new FileDialog(frame, "Open the coloring file", FileDialog.LOAD);
      while ( fd.getDirectory() == null && fd.getFile() == null) 
      {
        fd.show();
      }
      try {
        outputPdfDir = fd.getDirectory().toString();
        File got = new File( fd.getDirectory() + File.separator + fd.getFile());
        toRead = got.getAbsoluteFile().toURI().toURL().toString();
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
          shadings[p] = float(line[1]);
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
    keyPressed = false;
  }
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
  UpdateSimulation(false);
  keyPressed = eatKeyPress;
  if (renderThisFrame){
    endRecord();
  }
  textFont(txt);
  textSize(12);
  fill(0);
  textAlign(LEFT,TOP);
  text("Press any key to return.",0,0);
  
  if (keyPressed && (System.nanoTime()-gotHere)>.5e9){
    WHICH_SCREEN = 1;
    shadings = null;
  }
}





