
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import static java.lang.Math.sqrt;

// Ready for production use
public class PMFMe {

    public void calculateBoostPMF(String fileName) {
        final double c1 = 180.0 / Math.PI;
        String[] labels = {"D1", "D2", "X_phi", "X_psi", "D_phi", "D_psi", "F_phi", "F_psi", "D_chi1", "D_chi2", "F_chi1", "F_chi2"};
        String[] units = {"nm", "nm", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg", "deg"};
        double total = 0;
        double min = 0, max = 0;
        List<Double> valuesV = new ArrayList<>();
        List<List<Double>> valuesCV = new ArrayList<>();
        long N = 0;
        try (Scanner scan = new Scanner(new File(fileName))) {
            PrintWriter dihedrals = new PrintWriter("dihedrals.json");
            PrintWriter distances = new PrintWriter("distances.json");
            while (scan.hasNextLine()) {
                String line = scan.nextLine();
                if (line.indexOf("deltaV:") == 0) {
                    String[] nums = line.substring(line.indexOf("[")+1).replace("]", "").split(", ");
                    for (String v : nums) valuesV.add(Double.parseDouble(v));
                }
                if (line.indexOf("colvars:") == 0) {
                    String[] samples = line.substring(line.indexOf("(")+1).replace("]", "").replaceAll("\\(|\\)", " ").split(" ,  ");
                    for (String s : samples) {
                        String[] nums = s.split(", ");
                        List<Double> CV = new ArrayList<>();
                        for (String v : nums) {
                            CV.add(Double.parseDouble(v));
                        }
                        valuesCV.add(CV);
                    }
                    if (N++ % 20L == 0L) {
                        if (N == 1L) {
                            distances.print("[");
                            dihedrals.print("[");
                        }
                        else {
                            distances.print(",");
                            dihedrals.print(",");
                        }
                        distances.printf("[%f, %f]", valuesCV.get(valuesCV.size()-1).get(0), valuesCV.get(valuesCV.size()-1).get(1));
                        dihedrals.printf("[%f, %f, %f, %f, %f, %f, %f]", c1*valuesCV.get(valuesCV.size()-1).get(2), c1*valuesCV.get(valuesCV.size()-1).get(3),
                            c1*valuesCV.get(valuesCV.size()-1).get(4), c1*valuesCV.get(valuesCV.size()-1).get(5), c1*valuesCV.get(valuesCV.size()-1).get(6),
                            c1*valuesCV.get(valuesCV.size()-1).get(7), c1*valuesCV.get(valuesCV.size()-1).get(10));
                    }
                }
            }
            distances.print("]");
            dihedrals.print("]");
            distances.close();
            dihedrals.close();
            PrintWriter writer = new PrintWriter("output_dV.txt");
            min = valuesV.get(0); max = valuesV.get(0);
            for (Double d : valuesV) {
                total += d;
                if (d > max) max = d;
                if (d < min) min = d;
                writer.println(d);
            }
            writer.close();
        } catch (IOException e) {
            System.out.println("READ ERROR");
            return;
        }
        System.out.printf("Average boost potential %f kJ/mol (%f kcal/mol)\n", total/valuesV.size(), total/(4.184*valuesV.size()));
        System.out.printf("Range is %f to %f kJ/mol (%f to %f kcal/mol)\n", min, max, min/4.184, max/4.184);
        double mean = total/valuesV.size();
        double sum = 0.0;
        for (int i = 0; i < valuesV.size(); i++) {
            sum += Math.pow(valuesV.get(i) - mean, 2.0);
            valuesV.set(i, valuesV.get(i)-mean);
        }
        sum /= (valuesV.size() - 1);
        System.out.printf("The standard deviation is %f kJ/mol (%f kcal/mol)\n", sqrt(sum), sqrt(sum)/4.184);
        System.out.printf("%d points were incorporated\n", valuesV.size());
        for (int i = 0; i < 12; i++) {
            total = 0; sum = 0;
            max = valuesCV.get(0).get(i);
            min = valuesCV.get(0).get(i);
            for (int j = 0; j < valuesCV.size(); j++) {
                double value = (i > 1 ? 180.0*valuesCV.get(j).get(i)/Math.PI : valuesCV.get(j).get(i));
                if (value > max) max = value;
                if (value < min) min = value;
                total += value;
            }
            System.out.printf("\nAverage %s is %f %s\n", labels[i], total/valuesCV.size(), units[i]);
            System.out.printf("Range of %s is %f %s to %f %s\n", labels[i], min, units[i], max, units[i]);
            mean = total/valuesCV.size();
            ArrayList<Double> sorted = new ArrayList<>();
            for (int j = 0; j < valuesCV.size(); j++) {
                sum += (i > 1 ? Math.pow(180.0*valuesCV.get(j).get(i)/Math.PI - mean, 2.0) : Math.pow(valuesCV.get(j).get(i) - mean, 2.0));
                sorted.add((i > 1 ? 180*valuesCV.get(j).get(i)/Math.PI : valuesCV.get(j).get(i)));
                if (i > 1) valuesCV.get(j).set(i, 180.0*valuesCV.get(j).get(i)/Math.PI);
            }
            sorted.sort(Double::compare);
            sum /= (valuesCV.size() - 1);
            System.out.printf("The standard deviation of %s is %f %s\n", labels[i], sqrt(sum), units[i]);
            System.out.printf("The median of %s is %f %s\n", labels[i], sorted.get(sorted.size()/2), units[i]);
            System.out.printf("%d points were incorporated\n", valuesCV.size());
        }
        // compute histograms and obtain PMFs, 20 bins(?) for distances, 36 bins for dihedral angles
        // TO-DO: write distances at the end, determine ranges, and make bins for D1 and D2 (i = {0,1})
        //for (int i = 0; i < 2; i++) {

        //}
        for (int i = 2; i < 12; i++) {
            ArrayList<Double> bins = new ArrayList<>(36);
            ArrayList<Integer> hbins = new ArrayList<>(36);
            for (int j = 0; j < 36; j++) {
                bins.add(0.0);
                hbins.add(0);
            }
            for (int j = 0; j < valuesCV.size(); j++) {
                double value = (valuesCV.get(j).get(i) < 0 ? 360+valuesCV.get(j).get(i) : valuesCV.get(j).get(i));
                int bin = ((int)value)/10;
                bins.set(bin, bins.get(bin)+Math.exp(valuesV.get(j)/(4.184*0.616)));
                hbins.set(bin, hbins.get(bin)+1);
            }
            List<Double> pmf = new ArrayList<>();
            min = -Math.log(bins.get(0));
            for (int j = 0; j < 36; j++) {
                double value = -Math.log(bins.get(j));
                pmf.add(value);
                if (value < min) min = value;
            }
            System.out.printf("\n%s free energy histogram:\n", labels[i]);
            for (int j = 18; j < 36; j++) {
                System.out.printf("%s %f kT  (%d samples)\n", (j*10+5 - 360), pmf.get(j) - min, hbins.get(j));
            }
            for (int j = 0; j < 18; j++) {
                System.out.printf("%s %f kT  (%d samples)\n", (j*10+5), pmf.get(j) - min, hbins.get(j));
            }
        }
    }


    public static void main(String[] args) {
        if (args.length > 0) {
            PMFMe man = new PMFMe();
            man.calculateBoostPMF(args[0]);
        } else {
            System.out.println("Usage: java PMFMe fileName.out");
            System.out.println("fileName.out is an OpenMM output from CHARMM Python script with dV and colvar output lists");
        }
    }

}
