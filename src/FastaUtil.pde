
import javax.swing.JFileChooser;
import java.util.TreeSet;
import javax.swing.JOptionPane;
void TestFasta(int mode){

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
      toRead = got.getAbsoluteFile().toURI().toURL().toString();
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
        if (!buildString.matches("\\s+")){
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
      if (building!=null && !building.matches("\\s+")){
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
  int s = int((1+sqrt(1+8*Ncomp))/2); //AHAHAHHA fancy
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
void saveNodeStates(String outputFile){
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

    float PIX_TO_PD_FACTOR = (ZOOM/640.);

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

void GetCustomDistances(String customLoc){
  try {
    Scanner input = new Scanner(openStream(customLoc));
    int lineNum = 0;    
    while(input.hasNextLine()){
      input.nextLine(); 
      lineNum++;
    }
    input.close();
    input = new Scanner(openStream(customLoc));
    distanceMatrix = new float[lineNum][lineNum];
    for(int k = 0; k < lineNum; k++){
      String[] thisLine = input.nextLine().split("\\s+");
      for(int p = 0; p < lineNum; p++){
        distanceMatrix[k][p]=float(thisLine[p]);
      }
    }
    input.close();
  } 
  catch (Throwable e){
    e.printStackTrace();
  }
}
void pdScore2Mat(String[] lines){
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

void pdScoreMat(String[] lines){
  for(int k = 0; k < lines.length; k++){
    for(int p = 0; p < lines.length; p++){
      float val = PDScore(lines[p],lines[k]);
      System.out.printf("%.3f ",val);
    }
    System.out.println();
    System.err.println(k+" out of "+lines.length);
  }
}

void seqSimScore(String[] lines){
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
float seqSimDist(String a, String b){
  float diff = 0;
  for(int k = 0; k < a.length(); k++){
    if (k >= b.length() || a.charAt(k)!=b.charAt(k)){
        diff++;
    }
  }
  return diff / sqrt(a.length()*a.length()+b.length()*b.length()) * 25;
}

void pdScoreMultisearch(String[] lines){
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
      toRead = got.getAbsoluteFile().toURI().toURL().toString();
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
















