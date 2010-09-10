
float PDScore(String a, String b){
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
float PDScore0(String a, String b){
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
    0.354,	7.573,	11.294,	13.42,	-5.846,	6.599,	9.788,	9.655,	1.019,	-15.634,	-11.825,	10.762,	-10.585, -14.571, 7.662, 8.813,	3.012,	-13.11,	-6.245,	-12.135, 0,  }
  ,
  {
    3.762,	-10.135,	1.067,	-1.6,	4.885,	-5.166,	-7.861,	15.778,	-4.969,	1.993,	0.505,	-9.517,	-3.959,	-0.646,	8.029,	6.682,	4.127,	-5.222,	-1.6,	3.818,	0,  }
  ,
  {
    -11.036,	2.486,	2.718,	-0.325,	1.626,	-0.697,	-7.318,	-0.558,	0.953,	-2.045,	-6.157,	-1.022,	-3.601,	1.673,	9.456,	-0.348,	-0.348,	9.038,	9.874,	-4.345,	0,  }
  ,
  {
    -0.649,	-4.291,	1.963,	3.742,	9.397,	0.582,	2.611,	0.299,	4.657,	-3.243,	-4.557,	-5.405,	5.339,	-0.033,	-3.576,	-1.131,	-2.195,	1.38,	-1.597,	-3.26,	0,  }
  ,
  {
    2.828,	-5.687,	-0.859,	2.437,	-5.843,	-1.75,	4.734,	1.656,	-0.328,	-1.672,	3.219,	-0.422,	1.203,	3.25,	6,	-3.062,	-4.281,	4.64,	-1.422,	-4.672,	0,  }
};

float pdVal(char one, int att){
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
float PDScore2(String a, String b, int wSize){
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

