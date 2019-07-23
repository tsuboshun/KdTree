import java.io.*;
import infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov1;
import infodynamics.measures.continuous.kozachenko.EntropyCalculatorMultiVariateKozachenko;

class Test{
    public static void main(String args[]){
	int sample_N = 100000;
	String K = "100";
	
	MutualInfoCalculatorMultiVariateKraskov1 miCalc = null;
	EntropyCalculatorMultiVariateKozachenko enCalc = null;
	try{
	    miCalc = new MutualInfoCalculatorMultiVariateKraskov1();
	    miCalc.setProperty("k", K);
	    enCalc = new EntropyCalculatorMultiVariateKozachenko();
	}catch (InstantiationException e) {
	    e.printStackTrace();
	}catch (Exception e) {
	    e.printStackTrace();
	}

	double M_xy = 0;
	double H = 0;
	double xt, yt, zt;
	double[] xt_data = new double[sample_N];
	double[] yt_data = new double[sample_N];
	double[][] zt_data = new double[sample_N][1];
	FileReader fr = null;
	BufferedReader br = null;
	
	try{
	    miCalc.initialise();
	    miCalc.startAddObservations();
	    
	    fr = new FileReader("./test.txt");
	    br = new BufferedReader(fr);
	    String line;
	    String[] line_s;
	    int i = 0;
	    while ((line = br.readLine()) != null) {
		line_s = line.split("\\s");
		xt = Double.parseDouble(line_s[0]);
		yt = Double.parseDouble(line_s[1]);
		zt = Double.parseDouble(line_s[2]);
		xt_data[i] = xt;
		yt_data[i] = yt;
		zt_data[i][0] = zt;
		i++;
	    }

	    long start = System.currentTimeMillis();
	    enCalc.setObservations(zt_data);
	    H = enCalc.computeAverageLocalOfObservations();
	    long end = System.currentTimeMillis();
	    System.out.println("Shannon entropy " + H);
	    System.out.print(String.valueOf((end-start)/1000.0) + "sec\n");

	    start = System.currentTimeMillis();
	    miCalc.setObservations(xt_data, yt_data);
	    M_xy = miCalc.computeAverageLocalOfObservations();
	    end = System.currentTimeMillis();
	    System.out.println("Mutual information " + M_xy);
	    System.out.print(String.valueOf((end-start)/1000.0) + "sec\n");

	} catch (Exception e) {
	    e.printStackTrace();
	}	
    }
}
