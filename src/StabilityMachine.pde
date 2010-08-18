void getCoordsFile(float[][] coords, String file){
  String[] strings = loadStrings(file);
  for(int k = 0; k < strings.length; k++){
    String[] coord = strings[k].split("\\s+");
    if(coord.length < 5){
      continue;
    }
    coords[k][0] = new Float(coord[3]);
    coords[k][1] = new Float(coord[4]);
    float PIX_TO_PD_FACTOR = (ZOOM/640.);
    coords[k][0]/=PIX_TO_PD_FACTOR;
    coords[k][1]/=PIX_TO_PD_FACTOR;
  }
}
void rotateToBest(float[][] best, float[][] toModify, float[][] tmpMem, String outFile, PrintWriter dbg){
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
void rotateInto(float[][] one, float[][] dst, float theta){
  float ct = cos(theta); 
  float st = sin(theta);
  for(int k = 0; k < one.length; k++){
    dst[k][0] = ct*one[k][0]+st*one[k][1];
    dst[k][1] = -st*one[k][0]+ct*one[k][1];
  }
}
float calcBestRotation(float[][] one, float[][] two){
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
float sumDeviant(float[][] one, float[][] two){
  float toRet = 0;
  float dx, dy;
  for(int k = 0; k < one.length; k++){
    dx = one[k][0]-two[k][0];
    dy = one[k][1]-two[k][1];
    toRet += sqrt(dx*dx+dy*dy);
  }
  return toRet;
}




