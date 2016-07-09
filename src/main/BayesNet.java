/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package main;

import java.util.Random;

/**
 *
 * @author PeDeNRiQue
 */
public class BayesNet {
    
    static String[] hVistitAsia = {"A","P(A)"};
    static Double[][] mVisitAsia = {
        {1., 0.01},
        {0., 0.99}};
    
    static String[] hTuberculosis_VistiAsia = {"A","T","P(T|A)"};
    static Double[][] mTuberculosis_VisitAsia = {
        {1., 1., 0.05},
        {1., 0., 0.95},
        {0., 1., 0.01},
        {0., 0., 0.99}};
    
    static String[] hSmoke = {"S","P(S)"};
    static Double[][] mSmoke = {
        {1., 0.5},
        {0., 0.5}};
    
    static String[] hBronchitis = {"S","B","P(B|S)"};
    static Double[][] mBronchitis = {
        {1.,1.,0.6},
        {1.,0.,0.4},
        {0.,1.,0.3},
        {0.,0.,0.7}};
    
    static String[] hCancer = {"S","L","L|S"};
    static Double[][] mCancer = {
        {1.,1.,0.1},
        {1.,0.,0.9},
        {0.,1.,0.01},
        {0.,0.,0.99}};
    
    static String[] hTBorCancer = {"E","L","T","P(E|L,T"};
    static Double[][] mTBorCancer = {
        {1.,1.,1.,1.},
        {1.,1.,0.,1.},
        {1.,0.,1.,1.},
        {1.,0.,0.,0.},
        {0.,1.,1.,0.},
        {0.,1.,0.,0.},
        {0.,0.,1.,0.},
        {0.,0.,0.,1.}};
    
    static String[] hDispnea_Bronchitis_TB_Cancer = {"E","B","D","P(D|E,B)"};
    static Double[][] mDispnea_Bronchitis_TB_Cancer = {
        {1.,1.,1.,0.9},
        {1.,1.,0.,0.1},
        {1.,0.,1.,0.7},
        {1.,0.,0.,0.3},
        {0.,1.,1.,0.8},
        {0.,1.,0.,0.2},
        {0.,0.,1.,0.1},
        {0.,0.,0.,0.9}};
    
    static String[] hXRay = {"E","X","P(X|E)"};
    static Double[][] mXRay = {
        {1.,1.,0.98},
        {1.,0.,0.02},
        {0.,1.,0.05},
        {0.,0.,0.95}};
    
    final static int indexVisitAfrica = 0;
    final static int indexTuberculosis = 1;
    final static int indexTBorCancer = 2;
    final static int indexSmoker = 3;
    final static int indexCancer = 4;
    final static int indexBronchitis = 5;
    final static int indexDispnea = 6;
    final static int indexXRay = 7;
    final static int indexW = 8;
    final static int nIndexes = 9;
    
    static Boolean[] evidences = {
        true,//VisitAfrica
        true,//Tuberculosis
        true,//TBorCancer
        true,//Smoker
        true,//Cancer
        true,//Bronchitis
        true,//Dispnea
        true//XRay
    };
    
    public static void main(String[] args) {
        // TODO code application logic here
        int nSamples = 1000;
        Double[][] samples = new Double[nSamples][];
        
        for(int i = 0; i < nSamples; i++){
            //System.out.print(i+" | ");
            samples[i] = generateSample();
        }
        samples = joinSamples(samples);
        //smokerProbability(samples);
        //XRayProbability(samples);
        tuberculosisProbability(samples);
    }
    
    public static Double[][] joinSamples(Double[][] samples){
        int nSamples = samples.length;
        
        for(int i = 0; i < nSamples -1; i++){
            if(samples[i][0] == null){
                continue;
            }
            for(int j = i + 1; j < nSamples; j++){
                if(samples[j][0] == null){
                    continue;
                }
                if(isEquals(samples[i], samples[j])){
                    //System.out.println("IGUAIS: "+i+" "+j);
                    samples[i][indexW] +=  samples[j][indexW];
                    samples[j] = new Double[1];
                }                
            }
        }     
        
        for(int i = 0; i < nSamples; i++){
            showSample(samples[i]);
        }
        return samples;
    }
    
    public static Boolean isEquals(Double[] sample1,Double[] sample2){
    
        for(int c = 0; c < nIndexes - 1; c++){
            if(sample1[c].compareTo(sample2[c]) != 0){
                return false;
            }
        }
        
        return true;
    }
    
    public static void tuberculosisProbability(Double[][] samples){
        int nSamples = samples.length;
        Double total = 0.;
        Double tub = 0.;
        for(int i = 0; i < nSamples; i++){
            if(samples[i][0] == null){
                continue;
            }
            total += samples[i][indexW];
            if(samples[i][indexTuberculosis] == 0.0){
                tub += samples[i][indexW];
            }
        }
        
        System.out.printf("P(T) = %f\n",((1.0*tub)/total));
    }
    
    public static void XRayProbability(Double[][] samples){
        int nSamples = samples.length;
        Double total = 0.;
        Double x = 0.;
        for(int i = 0; i < nSamples; i++){
            if(samples[i][0] == null){
                continue;
            }
            total += samples[i][indexW];
            if(samples[i][indexXRay] == 0.0){
                x += samples[i][indexW];
            }
        }
        
        System.out.printf("P(X) = %f\n",(x / total));
    }
    
    public static void smokerProbability(Double[][] samples){
        int nSamples = samples.length;
        Double total = 0.;
        Double x = 0.;
        for(int i = 0; i < nSamples; i++){
            if(samples[i][0] == null){
                continue;
            }
            total += samples[i][indexW];
            if(samples[i][indexSmoker] == 0.0){
                x += samples[i][indexW];
            }
        }
        
        System.out.printf("P(X) = %f\n",(x / total));
    }
    
    public static Double[] generateSample(){
        Double[] sample = new Double[nIndexes];
        sample[indexW] = 1.; 
        sample = generateVisitAfrica(sample,false);
        sample = generateTuberculosis(sample,false);
        sample = generateSmoker(sample,false);
        sample = generateCancer(sample,false);
        sample = generateBronchiti(sample,false);
        sample = generateTBorCancer(sample,false);
        sample = generateXRay(sample,true);
        sample = generateDispnea(sample,true);
        //showSample(sample);
        
        return sample;
    }
    
    public static Double[] generateDispnea(Double[] sample,Boolean evidence){
        if(evidence == false){
            Double r = randomValue();
            if(sample[indexTBorCancer] == 1.0){
                if(sample[indexBronchitis] == 1.0){
                    if(r < mDispnea_Bronchitis_TB_Cancer[0][3]){
                        sample[indexDispnea] = 1.;
                    }else{
                        sample[indexDispnea] = 0.;
                    }
                }else{
                    if(r < mDispnea_Bronchitis_TB_Cancer[2][3]){
                        sample[indexDispnea] = 1.;
                    }else{
                        sample[indexDispnea] = 0.;
                    }
                }
            }else{
                if(sample[indexBronchitis] == 1.0){
                    if(r < mDispnea_Bronchitis_TB_Cancer[4][3]){
                        sample[indexDispnea] = 1.;
                    }else{
                        sample[indexDispnea] = 0.;
                    }
                }else{
                    if(r < mDispnea_Bronchitis_TB_Cancer[6][3]){
                        sample[indexDispnea] = 1.;
                    }else{
                        sample[indexDispnea] = 0.;
                    }
                }
            }
        }else{
            if(evidences[indexDispnea]){
                sample[indexDispnea] = 1.;
                if(sample[indexTBorCancer] == 1.0){
                    if(sample[indexBronchitis] == 1.0){
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[0][3];
                    }else{
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[2][3];
                    }
                }else{
                    if(sample[indexBronchitis] == 1.0){
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[4][3];
                    }else{
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[6][3];
                    }
                }
            }else{
                sample[indexDispnea] = 0.;
                if(sample[indexTBorCancer] == 1.0){
                    if(sample[indexBronchitis] == 1.0){
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[1][3];
                    }else{
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[3][3];
                    }
                }else{
                    if(sample[indexBronchitis] == 1.0){
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[5][3];
                    }else{
                        sample[indexW] *= mDispnea_Bronchitis_TB_Cancer[7][3];
                    }
                }
            }
        }        
        return sample;
    }
    
    public static Double[] generateXRay(Double[] sample,Boolean evidence){
        if(evidence == false){
            Double r = randomValue();
            if(sample[indexTBorCancer] == 1.0){
                if(r < mXRay[0][2]){
                    sample[indexXRay] = 1.;
                }else{
                    sample[indexXRay] = 0.;
                }
            }else{
                if(r < mXRay[2][2]){
                    sample[indexXRay] = 1.;
                }else{
                    sample[indexXRay] = 0.;
                }
            }
        }else{
            if(evidences[indexXRay]){
                sample[indexXRay] = 1.;
                if(sample[indexTBorCancer] == 1.0){
                    sample[indexW] *= mXRay[0][2];
                }else{
                    sample[indexW] *= mXRay[2][2];
                }
            }else{
                sample[indexXRay] = 0.;
                if(sample[indexTBorCancer] == 1.0){
                    sample[indexW] *= mXRay[1][2];
                }else{
                    sample[indexW] *= mXRay[3][2];
                }
            }
        }
        return sample;
    }
    
    public static Double[] generateTBorCancer(Double[] sample,Boolean evidence){
        if(evidence == false){
            if(sample[indexTuberculosis] == 1.0 || sample[indexCancer] == 1.0){
                sample[indexTBorCancer] = 1.;
            }else{
                sample[indexTBorCancer] = 0.;
            }
        }
        return sample;
    }
    
    public static Double[] generateBronchiti(Double[] sample,Boolean evidence){
        if(evidence ==  false){
            Double r = randomValue();
            if(sample[indexSmoker] == 1.0){
                if(r < mBronchitis[0][2]){
                    sample[indexBronchitis] = 1.;
                }else{
                    sample[indexBronchitis] = 0.;
                }
            }else{
                if(r < mBronchitis[2][2]){
                    sample[indexBronchitis] = 1.;
                }else{
                    sample[indexBronchitis] = 0.;
                }
            }
        }else{
            if(evidences[indexBronchitis]){
                sample[indexBronchitis] = 1.;
                if(sample[indexSmoker] == 1.0){
                    sample[indexW] *= mBronchitis[0][2];
                }else{
                    sample[indexW] *= mBronchitis[2][2];
                }
            }else{
                sample[indexCancer] = 0.;
                if(sample[indexSmoker] == 1.0){
                    sample[indexW] *= mBronchitis[1][2];
                }else{
                    sample[indexW] *= mBronchitis[3][2];
                }
            }
        }
        return sample;
    }
    
    public static Double[] generateCancer(Double[] sample,Boolean evidence){
        if(evidence ==  false){
            Double r = randomValue();
            if(sample[indexSmoker] == 1.0){
                if(r < mCancer[0][2]){
                    sample[indexCancer] = 1.;
                }else{
                    sample[indexCancer] = 0.;
                }
            }else{
                if(r < mCancer[2][2]){
                    sample[indexCancer] = 1.;
                }else{
                    sample[indexCancer] = 0.;
                }
            }
        }else{
            if(evidences[indexCancer]){
                sample[indexCancer] = 1.;
                if(sample[indexSmoker] == 1.0){
                    sample[indexW] *= mTuberculosis_VisitAsia[0][2];
                }else{
                    sample[indexW] *= mTuberculosis_VisitAsia[2][2];
                }
            }else{
                sample[indexCancer] = 0.;
                if(sample[indexSmoker] == 1.0){
                    sample[indexW] *= mTuberculosis_VisitAsia[1][2];
                }else{
                    sample[indexW] *= mTuberculosis_VisitAsia[3][2];
                }
            }
        }
        return sample;
    }
    
    public static Double[] generateSmoker(Double[] sample,Boolean evidence){
        if(evidence ==  false){
            Double r = randomValue();
            if(r < mSmoke[0][1]){
                sample[indexSmoker] = 1.;
            }else{
                sample[indexSmoker] = 0.;
            } 
        }else{
            if(evidences[indexSmoker]){
                sample[indexSmoker] = 1.;
                sample[indexW] *= mSmoke[0][1];
            }else{
                sample[indexSmoker] = 0.;
                sample[indexW] *= mSmoke[1][1];
            }
        }
        return sample;
    }
    
    public static Double[] generateTuberculosis(Double[] sample,Boolean evidence){
        if(evidence ==  false){
            Double r = randomValue();
            if(sample[indexVisitAfrica] == 1.0){
                if(r < mTuberculosis_VisitAsia[0][2]){
                    sample[indexTuberculosis] = 1.;
                }else{
                    sample[indexTuberculosis] = 0.;
                }
            }else{
                if(r < mTuberculosis_VisitAsia[2][2]){
                    sample[indexTuberculosis] = 1.;
                }else{
                    sample[indexTuberculosis] = 0.;
                }
            }
        }else{
            if(evidences[indexTuberculosis]){
                sample[indexTuberculosis] = 1.;
                if(sample[indexVisitAfrica] == 1.0){
                    sample[indexW] *= mTuberculosis_VisitAsia[0][2];
                }else{
                    sample[indexW] *= mTuberculosis_VisitAsia[2][2];
                }
            }else{
                sample[indexTuberculosis] = 0.;
                if(sample[indexVisitAfrica] == 1.0){
                    sample[indexW] *= mTuberculosis_VisitAsia[1][2];
                }else{
                    sample[indexW] *= mTuberculosis_VisitAsia[3][2];
                }
            }
        }
        return sample;
    }
    
    public static Double[] generateVisitAfrica(Double[] sample,Boolean evidence){
        if(evidence ==  false){
            Double r = randomValue();
            if(r < mVisitAsia[0][1]){
                sample[indexVisitAfrica] = 1.;
            }else{
                sample[indexVisitAfrica] = 0.;
            } 
        }else{
            if(evidences[indexVisitAfrica]){
                sample[indexVisitAfrica] = 1.;
                sample[indexW] *= mVisitAsia[0][1];
            }else{
                sample[indexVisitAfrica] = 0.;
                sample[indexW] *= mVisitAsia[1][1];
            }
        }
        return sample;
    }
    
    public static Double randomValue(){
        Double rangeMin = 0.;
        Double rangeMax = 1.;
        
        Random r = new Random();
        double randomValue = rangeMin + (rangeMax - rangeMin) * r.nextDouble();
        //System.out.println(randomValue+"");
        return randomValue;
    }
    
    public static void showSample(Double[] sample){
        if(sample[0] == null){
            return;
        }
        for(int i = 0; i < nIndexes; i++){
            if(i < 8){
                System.out.print(sample[i]+"  ");
            }else{
                System.out.printf(" | %f\n",sample[i]);
            }
        }
    }
}
