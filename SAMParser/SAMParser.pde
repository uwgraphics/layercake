/*public enum SAMLine{
 QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL
 }*/
final static int QNAME = 0;
final static int FLAG = 1;
final static int RNAME = 2;
final static int POS = 3;
final static int MAPQ = 4;
final static int CIGAR = 5;
final static int RNEXT = 6;
final static int PNEXT = 7;
final static int TLEN = 8;
final static int SEQ = 9;
final static int QUAL = 10;

final static int A = 0;
final static int T = 1;
final static int C = 2;
final static int G = 3;

final static char[] bps = {'A','T','C','G'};
int[][] seq;
String name;
boolean printed = false;


//parsing info
int minrt = 999999;
int maxrt = 0;
int rows = 0;
HashMap names = new HashMap();
float[][][] variants;
double[][][] pvalues;
int[][] coverage;
char[][] refgen;
boolean loaded = false;

void setup() {
  size(400,100);
  background(255);
  fill(0);
  stroke(0);

}

void draw(){
  if(!loaded){
    selectFolder("Select a folder of SAM files:", "parseSAMFolder");
    //parseSAMFolder(new File(dataPath("")));
    loaded = true;  
  }
  
}

boolean isHeaderLine(String header) {
  if (header.equals("@HD") || header.equals("@SQ")||header.equals("@RG")||header.equals("@PG"))
    return true;
  return false;
}

void readSAMLine(String[] tokens) {

  if (!printed) {
    System.out.println(tokens.length+" fields present in this line:");
    for (String token:tokens) {
      System.out.println(token);
    }
    printed = true;
  }
}

int checkRowsSAM(String filename,int curdex) {
  System.out.println("Checking "+filename+" for individual sequences");
  BufferedReader reader = createReader(filename);
  String aline = "";
  String[] tokens;
  int index = curdex;
  while (aline!=null) {
    try {
      aline = reader.readLine();
      if (aline!=null) {
        tokens = splitTokens(aline, "\t");

        if (tokens!=null && tokens.length>0 && !isHeaderLine(tokens[0])) {
          if (!names.containsKey(filename+"_"+tokens[2])) {
            System.out.println("Found sequence "+tokens[2]);
            names.put(filename+"_"+tokens[2],index);
            index++;
          }
          if (int(tokens[3])<minrt) {
            minrt = int(tokens[3]);
          }
          if (int(tokens[3])+tokens[9].length()>maxrt) {
            maxrt = int(tokens[3])+tokens[9].length();
          }
        }//end if
        
      }
    }//end try
    catch(IOException e) {
      e.printStackTrace();
    }//end catch
  }//end while
  System.out.println("Finished checking "+filename+", found "+(index-curdex)+" sequences");
  return index;
}//end checkRowsSAM


void parseSAMFolder(File folder) {
  if(folder==null){
    exit();
  }
  else{
  ArrayList<String> samFiles = new ArrayList<String>(); 
  String ext;
  rows = 0;
  if (folder!=null) {
      for(File file: folder.listFiles()){
         ext = file.getPath().substring(file.getPath().lastIndexOf('.')+1);
         if(ext.equals("sam")){
           samFiles.add(file.getPath());
         }
      }
      
    }
    
    float dx = 350.0/samFiles.size();
    text("Found "+samFiles.size()+" SAM files:",5,10);
    text("Pre-Processing:",5,20);
    fill(255);
    stroke(0);
    rect(25,40,350,20);
    float x = 25;
        
    for(String aFile:samFiles){
      rows = checkRowsSAM(aFile,rows);
      noStroke();
      fill(0,0,255);
      rect(x,40,dx,20);
      x+=dx;
      
    }
    

    int cols = maxrt-minrt;
    System.out.println("Making "+cols+" bins");
    variants = new float[rows][cols][4];
    pvalues = new double[rows][cols][4];
    refgen = new char[rows+1][cols];
    coverage = new int[rows][cols];
    
    for(int i = 0;i<variants.length;i++){
      for (int j = 0;j<variants[i].length;j++) {
        variants[i][j] = new float[] {
        0, 0, 0, 0
        };
        pvalues[i][j] = new double[] {
          0, 0, 0, 0
        };
        coverage[i][j] = 0;
      }//end for
    }
    
    System.out.println("Loading data");
    fill(0);
    text("Loading:",5,30);
    
    fill(255);
    stroke(0);
    rect(25,65,350,20);
    
    int numDone = 0;
    for(String aFile:samFiles){
      loadOneSAM(aFile);
      noStroke();
      fill(0,0,255);
      rect(25+(dx*numDone),65,dx,20);
      numDone++;
    }
    
    
    float vsum = 0;
    float curmax = 0;
    char curref = 'N';
    for(int i = 0;i<rows;i++){
      for(int j = 0;j<variants[i].length;j++){
        vsum = sum(variants[i][j]);
        curref = 'N';
        curmax = 0;
        for (int k=0;k<variants[i][j].length;k++){
          if(variants[i][j][k]>curmax){
           curmax = variants[i][j][k];
            curref = bps[k];
          }
          if(variants[i][j][k]>0){
            pvalues[i][j][k] /= variants[i][j][k];
          } 
          variants[i][j][k] /= coverage[i][j];

        }
        refgen[i+1][j] = curref; 
      }
    }
    fill(0);
    text("Writing...",5,95);
    System.out.println("Writing out file");
    int numBins = 0;
    for(int i = 0;i<coverage.length;i++){
      for(int j= 0;j<coverage[i].length;j++){
        if(coverage[i][j]>0){
          numBins++;
        }
      }
    }
    
    System.out.println("Should have "+numBins+" non empty bins");
  //  prepMatrix(samFiles.get(0));
    makeCSV(folder);
    exit();
  }
}


String[] toMetaNames(){
  Object[] temp = names.keySet().toArray();
  String fullName;
  String[] fileNames = new String[temp.length];
  String[] seqNames = new String[temp.length];
  String[] metadata = new String[temp.length];

  int delimiter;
  for(int i = 0;i<metadata.length;i++){
    int index = ((Integer)(names.get(temp[i]))).intValue();
    fullName = (String)temp[i];
    delimiter = fullName.indexOf(".sam_");
    fileNames[index] = fullName.substring(0,delimiter);
    seqNames[index] = fullName.substring(delimiter+5);
  }
  
  for(int i = 0;i<metadata.length;i++){
    if(numFound(seqNames,seqNames[i])>1){
      if(fileNames[i].lastIndexOf("/")!=-1){
        metadata[i] = fileNames[i].substring(fileNames[i].lastIndexOf("/")+1);
      }
      else if(fileNames[i].lastIndexOf("\\")!=-1){
        metadata[i] = fileNames[i].substring(fileNames[i].lastIndexOf("\\")+1);
      }
      else{
        metadata[i] = fileNames[i];
      }
    }
    else{
      metadata[i] = seqNames[i];
    }
  }
  return metadata;
  
}

int numFound(String[] source ,String target){
  int found = 0;
  for(String aString:source){
    if(aString.equals(target)){
      found++;
    }
  }
  return found;
}


void makeCSV(File folder) {
  String[] metadata = toMetaNames();
  
  System.out.println("Making file with "+metadata.length+" sequences, maximum of "+variants[0].length+" bps per sequence");
  PrintWriter file = createWriter(folder.getAbsolutePath()+"/LayerCakeInput.csv");
  file.println("Sequence Name,Reference Nucleotide(s),Variant Nucleotide(s),Nucleotide,Coverage,Variant Frequency,Variant P-Value (approximate)"); 
    for(int i = 0;i<variants.length;i++){
      for(int j = 0;j<variants[i].length;j++){
        if(coverage[i][j]>0){
          for(int k=0;k<variants[i][j].length;k++){
         //   if(variants[i][j][k]>0)
              file.println(metadata[i]+","+refgen[i+1][j]+","+bps[k]+","+(j+minrt)+","+coverage[i][j]+","+floatToPercent(variants[i][j][k])+","+pvalues[i][j][k]);
          }
        }
      }
    }
  file.flush();
  file.close();
  System.out.println("Done Writing");
}


String floatToPercent(float num){
  return (round(num*10000)/100.0)+"%";
}

float sum(float[] anArray){
  float val = 0;
  for(float aVal:anArray){
    val+=aVal;
  }
  return val;
}

float sum(int[] anArray){
  float val = 0;
  for(int aVal:anArray){
    val+=aVal;
  }
  return val;
}

boolean loadOneSAM(String filename) {
  // System.err.println("Loading "+filename);
  int index = 0;

  BufferedReader reader = createReader(filename);
  String aline = "";
  String[] tokens;
  int c;
  char[] reads;
  boolean inBounds;
  int extent;
  int bp_i;
  double phred;
  int numBins = 0;
  while (aline!=null) {
    try {
      aline = reader.readLine();
      if (aline!=null) {
        tokens = splitTokens(aline, "\t");
        if (tokens!=null && tokens.length>0 && !isHeaderLine(tokens[0])) {
          index = ((Integer)(names.get(filename+"_"+tokens[2]))).intValue();
          c = int(tokens[3])-minrt;
          extent = tokens[9].length();
          if (c+extent >= 0) {
            reads = tokens[9].toCharArray();
            for (int i = 0;i<reads.length;i++) {
              inBounds = (c+i>=0 && c+i<variants[index].length);
              if (inBounds) {
                bp_i = bpToIndex(reads[i]);

                if (!(tokens[10].charAt(i)=='*')) {

                    

                  coverage[index][c+i]++;
                  phred = int(tokens[10].charAt(i))-33d;
                  if(bp_i>=0){
                    pvalues[index][c+i][bp_i]+= Math.pow(10d,phred/-10);
                    variants[index][c+i][bp_i]+=1.0;
                  }
                  else{
                    System.err.print(reads[i]);
                  }
                  numBins++;
                }//end if
              }//end if
            }//end for
          }//end if
        }//end if
      }//end if
    }//end try
    catch(IOException e) {
      e.printStackTrace();
    }//end catch
  }//end while

  System.out.println("Successfully loaded "+filename+" with "+numBins+" reads");
  return true;
}//end loadOneSAM

int bpToIndex(char bp) {
  int index = -1;
  switch(bp) {
  case 'A':
  case 'a':
    index = 0;
    break;

  case 'T':
  case 't':
    index =  1;
    break;

  case 'C':
  case 'c':
    index =  2;
    break;

  case 'G':
  case 'g':
    index =  3;
    break;

  default:
    break;
  }//end switch
  return index;
}


void prepMatrix(String filename){
  System.err.println("Prepping Matrix");
  float[][] cooccurance = new float[maxrt-minrt][maxrt-minrt];
    float[][] coverages = new float[maxrt-minrt][maxrt-minrt];
  for(int i = 0;i<cooccurance.length;i++){
    for(int j=0;j<cooccurance[i].length;j++){
      cooccurance[i][j] = 0;
      coverages[i][j] = 0;
    }
  }
    int index = 0;

  BufferedReader reader = createReader(filename);
  String aline = "";
  String[] tokens;
  int c;
  char[] reads;
  boolean inBounds;
  int extent;
  int bp_i;
  boolean[] isVariant;
  
  
  while (aline!=null) {
    try {
      aline = reader.readLine();
      if (aline!=null) {
        tokens = splitTokens(aline, "\t");
        if (tokens!=null && tokens.length>0 && !isHeaderLine(tokens[0])) {
          index = ((Integer)(names.get(filename+"_"+tokens[2]))).intValue();
          c = int(tokens[3])-minrt;
          extent = tokens[9].length();
          if (c+extent >= 0) {
            reads = tokens[9].toCharArray();
            isVariant = new boolean[reads.length];
            for(int i = 0;i<isVariant.length;i++){
              isVariant[i]= false;
            }
            for (int i = 0;i<reads.length;i++) {
              inBounds = (c+i>=0 && c+i<variants[index].length);
              if (inBounds) {
                bp_i = bpToIndex(reads[i]);
                isVariant[i] = ((refgen[index+1][c+i]!=reads[i])&&(tokens[10].charAt(i)!='*'));
              }//end if
            }//end for
            for (int i = 0;i<reads.length;i++) {
              inBounds = (c+i>=0 && c+i<variants[index].length);
              if (isVariant[i] && inBounds && !(tokens[10].charAt(i)=='*')) {
                bp_i = bpToIndex(reads[i]);
                   for(int j = 0;j<reads.length;j++){
                     if(isVariant[j]){
                       cooccurance[c+i][c+j]++;
                     }//end if
                   }//end for
              }//end if
            }//end for
          }//end if
        }//end if
      }//end if
    }//end try
    catch(IOException e) {
      e.printStackTrace();
    }//end catch
  }//end while
  
  System.err.println("Drawing Matrix");
  size(cooccurance.length,cooccurance.length);
  noStroke();
  float divisor;
  for(int i = 0;i<cooccurance.length;i++){
    for(int j=0;j<cooccurance[i].length;j++){
      if(bpToIndex(refgen[1][i])>=0)
       divisor = (1.0-variants[0][i][bpToIndex(refgen[1][i])])*coverage[0][i];
      else
       divisor = coverage[0][i]; 
      if(divisor>0){
        cooccurance[i][j]/= (divisor);
      }
      
      
      fill(255*(cooccurance[i][j]),0,0);
      
      rect(j,i,1,1);
    }
  }
  save(dataPath("CoOccuranceMatrix.png"));
  System.err.println("Printing out matrix");

  
  PrintWriter file = createWriter(dataPath("CoOccurranceMatrix.csv"));
  for(int i = 0;i<10;i++){//cooccurance.length;i++){
    for(int j = 0;j<cooccurance[i].length-1;j++){
      file.print(cooccurance[i][j]+",");
    }
    file.print(cooccurance[i][cooccurance[i].length-1]+"");
    file.println("");
  }
  file.flush();
  file.close();
  

  
  System.err.println("Done");
  
}

