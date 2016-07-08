/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package entities;

import java.util.List;

/**
 *
 * @author PeDeNRiQue
 */
public class Table {
    private String[] headers;
    private Double[][] matrix;

    public Table(String[] headers,Double[][] matrix){
        this.headers = headers;
        this.matrix = matrix;
    }

    public String[] getHeaders() {
        return headers;
    }

    public void setHeaders(String[] headers) {
        this.headers = headers;
    }
    
    public Double[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(Double[][] matrix) {
        this.matrix = matrix;
    }
    
    public int size(){
        return matrix[0].length;
    }
    
    public int hasHeader(String header){
        for(int i = 0; i < this.headers.length; i++){
            if(this.headers[i].equals(header)){
                return i;
            }
        }
        return -1;
    }
}
