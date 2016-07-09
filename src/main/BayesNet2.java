/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package main;

/**
 *
 * @author PeDeNRiQue
 */
public class BayesNet2 {
    static String[] hVistitAsia = {"A","P(A)"};
    static Double[][] mVisitAsia = {{1., 0.01},
                                    {0., 0.99}};
    
    static String[] hTuberculosis_VistitAsia = {"A","T","P(T|A)"};
    static Double[][] mTuberculosis_VisitAsia = {{1., 1., 0.05},
                                          {1., 0., 0.95},
                                          {0., 1., 0.01},
                                          {0., 0., 0.99}};
    
    static String[] hSmoke = {"S","P(S)"};
    static Double[][] mSmoke = {{1., 0.5},
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
        {1.,0.,0.,0.}};
    
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
    
    public static void main(String[] args) {
        
        Double[][] result = calculateLBS(mSmoke,mBronchitis,mCancer);
        result = eliminateS_resultLB(result);
        result = calculateELTB(mTBorCancer,result);
        result = normalizeELTB(result); 
        result = eliminateL_resultETB(result);
        result = calculateDEBT(mDispnea_Bronchitis_TB_Cancer,result);
        result = eliminateB_resultEDT(result);
        result = calculatEDTX(mXRay,result);
        result = eliminateE_resultDTX(result);
        result = eliminateX_resultDT(result);
        
        result = calculateADT(calculateTA(mVisitAsia,mTuberculosis_VisitAsia),result);
        result = eliminateAD_resultT(result);
    }
    
    public static Double[][] eliminateAD_resultT(Double[][] matrix_DTX){
        System.out.println("\nELIMINANDO X -> P(D,T)\n");
        Double[][] matrix_T = new Double[2][2];
        
        matrix_T[0][1] = 0.;
        matrix_T[1][1] = 0.;
        
        for(int i = 0; i < 8; i++){
            matrix_T[i % 2][1] +=  matrix_DTX[i][3];            
        }
        
        System.out.println("T  | P(T)");
        System.out.println("+t | "+matrix_T[0][1]);
        System.out.println("-t | "+matrix_T[1][1]);
        return matrix_T;
    }
    
    public static Double[][] calculateADT(Double[][] matrix_AT,Double[][] matrix_DT){
        System.out.println("\nCALCULANDO P(A,D,T)\n");
        Double[][] matrix_ADT = new Double[8][4];
        int i1,i2;
        
        for(int i = 0; i < 8; i++){
            
            if(i < 4){
                matrix_ADT[i][0] = 1.;
            }else{
                matrix_ADT[i][0] = 0.;
            }
            
            if(i == 0 || i == 1 || i == 4 || i == 5){
                matrix_ADT[i][1] = 1.;
            }else{
                matrix_ADT[i][1] = 0.;
            }
            
            if(i % 2 == 0){
                matrix_ADT[i][2] = 1.;
            }else{
                matrix_ADT[i][2] = 0.;
            }
            
            if(i < 4){
                if(i % 2 == 0){
                    i1 = 0;
                }else{
                    i1 = 1;
                }
            }else{
                if(i % 2 == 0){
                    i1 = 2;
                }else{
                    i1 = 3;
                }
            }
            
            if(i < 4){
                i2 = i;
            }else{
                i2 = i - 4;
            }
            
            matrix_ADT[i][3] = matrix_AT[i1][2] *  matrix_DT[i2][2];
        }
        Double total = 0.;
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 4; j++){
                if(j == 3)
                    total += matrix_ADT[i][j];
                //System.out.printf("%f ",matrix_ADT[i][j]);
            }
            //System.out.println("");
        }
        System.out.println("\nNORMALIZANDO\n");
        Double total2 = 0.;
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 4; j++){
                if(j == 3){
                    matrix_ADT[i][j] = matrix_ADT[i][j] / total;
                    total2 += matrix_ADT[i][j];
                }
                System.out.printf("%f ",matrix_ADT[i][j]);
            }
            System.out.println("");
        }
        //System.out.println("TOTAL: "+total2);
        
        return matrix_ADT;
    }
    
    public static Double[][] eliminateX_resultDT(Double[][] matrix_DTX){
        System.out.println("\nELIMINANDO X -> P(D,T)\n");
        Double[][] matrix_DT = new Double[4][3];
    
        for(int i = 0; i < 4; i++){
            if(i < 2){
                matrix_DT[i][0] = 1.;
            }else{
                matrix_DT[i][0] = 0.;
            }
            if(i % 2 == 0){
                matrix_DT[i][1] = 1.;
            }else{
                matrix_DT[i][1] = 0.;
            }
            
            matrix_DT[i][2] = matrix_DTX[i*2][3] + matrix_DTX[(i*2)+1][3];
        }
        
        for(int i = 0; i < 4 ;i++){
            for(int j = 0; j < 3; j++){
                System.out.print(matrix_DT[i][j]+" ");
            }
            System.out.println("");
        }
        
        return matrix_DT;
    }
    

    public static Double[][] eliminateE_resultDTX(Double[][] matrix_EDTX){
        System.out.println("\nELIMINANDO E -> P(D,T,X)\n");
        Double[][] matrix_DTX = new Double[8][4];
     
        for(int i = 0; i < 8; i++){
            for(int j = 0; j < 4; j++){
                matrix_DTX[i][j] = matrix_EDTX[i][j+1];
            }
        }
        
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 4; j++){
                System.out.printf("%f ",matrix_DTX[i][j]);
            }
            System.out.println("");
        }
        
        return matrix_DTX;
    }
    
    public static Double[][] calculatEDTX(Double[][] matrix_XE,Double[][] matrix_EDT){
        System.out.println("\nCALCULANDO P(E,D,T,X)\n");
        Double[][] matrix_EDTX = new Double[8][5];
        int i1,i2;
        
        for(int i = 0; i < 8; i++){
            matrix_EDTX[i][0] = 1.;
            
            if(i < 4){
                matrix_EDTX[i][1] = 1.;
            }else{
                matrix_EDTX[i][1] = 0.;
            }
            
            if(i == 0 || i == 1 || i == 4 || i == 5){
                matrix_EDTX[i][2] = 1.;
            }else{
                matrix_EDTX[i][2] = 0.;
            }
            
            if(i % 2 == 0){
                matrix_EDTX[i][3] = 1.;
            }else{
                matrix_EDTX[i][3] = 0.;
            }
            
            if(i % 2 == 0){
                i1 = 0;
            }else{
                i1 = 1;
            }
            i2 = i / 2;
            
            matrix_EDTX[i][4] = matrix_XE[i1][2] *  matrix_EDT[i2][3];
        }
        
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 5; j++){
                System.out.print(matrix_EDTX[i][j]+" ");
            }
            System.out.println("");
        }
        return matrix_EDTX;
    }
    
    public static Double[][] eliminateB_resultEDT(Double[][] matrix_EBDT){
        System.out.println("\nELIMINANDO B -> P(E,D,T)\n");
        Double[][] matrix_EDT = new Double[4][4];
        
        for(int i = 0; i < 4; i++){
            matrix_EDT[i][0] = 1.;
            
            if(i < 2){
                 matrix_EDT[i][1] = 1.;
            }else{
                 matrix_EDT[i][1] = 0.;
            }
            
            if(i % 2 == 0){
                 matrix_EDT[i][2] = 1.;
            }else{
                 matrix_EDT[i][2] = 0.;
            }
            
            matrix_EDT[i][3] = matrix_EBDT[2*i][4] + matrix_EBDT[(2 * i)+1][4];
        }
        for(int i = 0; i < 4 ;i++){
            for(int j = 0; j < 4; j++){
                System.out.print(matrix_EDT[i][j]+" ");
            }
            System.out.println("");
        }
        return matrix_EDT;
    }
    
    public static Double[][] calculateDEBT(Double[][] matrix_DEB,Double[][] matrix_ETB){
        System.out.println("\nCALCULANDO DEBT\n");
        Double[][] matrix_DEBT = new Double[8][5];
        int i1,i2;
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                matrix_DEBT[i][1] = 1.;
            }else{
                matrix_DEBT[i][1] = 0.;
            }
            
            if(i == 0 || i == 1 || i == 4 || i == 5){
                matrix_DEBT[i][2] = 1.;
            }else{
                matrix_DEBT[i][2] = 0.;
            }
            
            if(i % 2 == 0){
                matrix_DEBT[i][3] = 1.;
            }else{
                matrix_DEBT[i][3] = 0.;
            }
            
            matrix_DEBT[i][0] = 1.;
            
            
            if(i < 4){
                if(i % 2 == 0){
                    i1 = 0;
                }else{
                    i1 = 2;
                }
            }else{
                if(i % 2 == 0){
                    i1 = 1;
                }else{
                    i1 = 3;
                }
            }
            
            if(i == 0|| i == 1 || i == 4 || i == 5){
                if(i % 2 == 0){
                    i2 = 0;
                }else{
                    i2 = 1;
                }
            }else{
                if(i % 2 == 0){
                    i2 = 2;
                }else{
                    i2 = 3;
                }
            }
            
            matrix_DEBT[i][4] = matrix_DEB[i1][3] * matrix_ETB[i2][3];
            //System.out.println(matrix_DEB[i1][3]+" "+matrix_ETB[i2][3]+" = "+matrix_DEBT[i][4]);
        }
        
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 5; j++){
                
                System.out.print(matrix_DEBT[i][j]+" ");
            }
            System.out.println("");
        }
        
        return matrix_DEBT;
    }
    
    public static Double[][] eliminateL_resultETB(Double[][] matrix_LETB){
        System.out.println("\nELIMINANDO L -> P(E,T,B)\n");
        Double[][] matrix_ETB = new Double[4][4];
        
        for(int i = 0; i < 4; i++){
            matrix_ETB[i][0] = 1.;
            
            if(i == 0 || i == 1 || i == 4 || i == 5){
                matrix_ETB[i][1] = 1.;
            }else{
                matrix_ETB[i][1] = 0.;
            }
            
            if(i % 2 == 0){
                matrix_ETB[i][2] = 1.;
            }else{
                matrix_ETB[i][2] = 0.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            matrix_ETB[i][3] = matrix_LETB[i][4] + matrix_LETB[i+4][4];
        }
        for(int i = 0; i < 4 ;i++){
            for(int j = 0; j < 4; j++){
                System.out.print(matrix_ETB[i][j]+" ");
            }
            System.out.println("");
        }
        return matrix_ETB;
    }
    
    public static Double[][] normalizeELTB(Double[][] matrix_ELT){
        System.out.println("\nNORMALIZANDO VALORES DE P(E,L,T)\n");
        Double total = 0.;
        
        for(int i = 0; i < 8; i++){
            total = total + matrix_ELT[i][4];
        }
        
        for(int i = 0; i < 8; i++){
            matrix_ELT[i][4] = matrix_ELT[i][4] / total;
        }
        
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 5; j++){
                System.out.print(matrix_ELT[i][j]+" ");
            }
            System.out.println("");
        }
        
        return matrix_ELT;
    }
    
    public static Double[][] calculateELTB(Double[][] matrix_ELT,Double[][] matrix_LB){
        System.out.println("\nCALCULANDO P(E,L,T,B)\n");
        Double[][] matrix_ELTB = new Double[8][5];
        int i1,i2;
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                matrix_ELTB[i][1] = 1.;
            }else{
                matrix_ELTB[i][1] = 0.;
            }
            
            if(i == 0 || i == 1 || i == 4 || i == 5){
                matrix_ELTB[i][2] = 1.;
            }else{
                matrix_ELTB[i][2] = 0.;
            }
            
            if(i % 2 == 0){
                matrix_ELTB[i][3] = 1.;
            }else{
                matrix_ELTB[i][3] = 0.;
            }
            
            matrix_ELTB[i][0] = 1.;
            
            i1 = i / 2;
            
            if(i < 4){
                if(i % 2 == 0){
                    i2 = 0;
                }else{
                    i2 = 1;
                }
            }else{
                if(i % 2 == 0){
                    i2 = 2;
                }else{
                    i2 = 3;
                }
            }
            
            matrix_ELTB[i][4] = matrix_ELT[i1][3] * matrix_LB[i2][2];
        }
        
        for(int i = 0; i < 8 ;i++){
            for(int j = 0; j < 5; j++){
                System.out.print(matrix_ELTB[i][j]+" ");
            }
            System.out.println("");
        }
        
        return matrix_ELTB;
    }
    
    public static Double[][] eliminateS_resultLB(Double[][] matrix_LBS){
        Double[][] matrix_LB = new Double[4][3];
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                matrix_LB[i][0] = 1.;
            }else{
                matrix_LB[i][0] = 0.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            if(i % 2 == 1){
                matrix_LB[i][1] = 0.;
            }else{
                matrix_LB[i][1] = 1.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            matrix_LB[i][2] = matrix_LBS[i*2][3] + matrix_LBS[(i*2)+1][3];
        }
        
        
        for(int i = 0; i < 4 ;i++){
            for(int j = 0; j < 3; j++){
                System.out.print(matrix_LB[i][j]+" ");
            }
            System.out.println("");
        }
        
        return matrix_LB;
    }
     
    public static Double[][] calculateLBS(Double[][] matrix_S, Double[][] matrix_S_B,Double[][] matrix_S_L){
        Double[][] matrix_LBS = new Double[8][4];
        
        for(int i = 0; i < 8; i++){
            if(i < 4){
                matrix_LBS[i][0] = 1.;
            }else{
                matrix_LBS[i][0] = 0.;
            }
        }
        
        for(int i = 0; i < 8; i = i + 2){
            if((i/2)  % 2 == 0){
                matrix_LBS[i][1] = 1.;
                matrix_LBS[i+1][1] = 1.;
            }else{
                matrix_LBS[i][1] = 0.;
                matrix_LBS[i+1][1] = 0.;
            }
        }
        
        for(int i = 0; i < 8; i++){
            if(i % 2 == 0){
                matrix_LBS[i][2] = 1.;
            }else{
                matrix_LBS[i][2] = 0.;
            }
        }
        
                
        int i1,i2,i3;
        
        for(int i = 0; i < 8; i++){
            if(i % 2 == 0){
                i1 = 0;
            }else{
                i1 = 1;
            }
            
            if(i1 == 0){
                if(i < 4){
                    i2 = 0;
                }else{
                    i2 = 1;
                }
            }else{
                if(i < 4){
                    i2 = 2;
                }else{
                    i2 = 3;
                }
            }
            
            if(i1 == 0){
                if((i/2)  % 2 == 0){
                    i3 = 0;
                }else{
                    i3 = 1;
                }
            }else{
                if((i/2)  % 2 == 0){
                    i3 = 2;
                }else{
                    i3 = 3;
                }
            }
            
            matrix_LBS[i][3] = matrix_S[i1][1] * matrix_S_B[i3][2] * matrix_S_L[i2][2];
        }
        
//        for(int i = 0; i < 8; i++){
//            for(int j = 0; j < 4; j++){
//                
//                System.out.print(matrix_LBS[i][j]+" ");
//            }
//            System.out.println("");
//        }
        
        return matrix_LBS;
    }
    
    public static Double[][] calculateTA(Double[][] matrix_A, Double[][] matrix_T_A){
        System.out.println("\nCALCULANDO TA\n");
        Double[][] matrix_TA = new Double[4][3];
        int index1,index2;
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                matrix_TA[i][0] = 1.;
            }else{
                matrix_TA[i][0] = 0.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            if(i % 2 == 1){
                matrix_TA[i][1] = 0.;
            }else{
                matrix_TA[i][1] = 1.;
            }
        }
        
        for(int i = 0; i < 4; i++){
            if(i < 2){
                    index1 = 0;
            }else{
                index1 = 1;
            }
            matrix_TA[i][2] = matrix_A[index1][1] * matrix_T_A[i][2]; 
        }
               
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 3; j++){
                
                System.out.print(matrix_TA[i][j]+" ");
            }
            System.out.println("");
        }
        
        return matrix_TA;
    }
}
