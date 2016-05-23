package dictionary;

import Jama.Matrix;

public class MatrixForDic extends Matrix {

	public MatrixForDic(double[][] arg0) {
		super(arg0);
	}
	public MatrixForDic(int arg0, int arg1) {
		super(arg0, arg1);
	}
	public MatrixForDic(Matrix matrix) {
		super(matrix.getArrayCopy());
	}

	public double[] getRowCopy(int index) {
		if((index < 0) || (getRowDimension() <= index))
			return null;
		double[] row = new double[getColumnDimension()];
		for(int i = 0; i < row.length; i++)
			row[i] = this.get(index, i);
		return row;
	}
	public double[] getColumnCopy(int index) {
		if((index < 0) || (getColumnDimension() <= index))
			return null;
		double[] column = new double[getRowDimension()];
		for(int i = 0; i < column.length; i++)
			column[i] = this.get(i, index);
		return column;
	}
	public MatrixForDic setRow(int index, double[] row) {
		if((index < 0) || (getRowDimension() <= index) ||
				(row.length != getColumnDimension()))
			return null;
		for(int i = 0; i < row.length; i++)
			this.set(index, i, row[i]);
		return this;
	}
	public MatrixForDic setColumn(int index, double[] column) {
		if((index < 0) || (getColumnDimension() <= index) ||
				(column.length != getRowDimension()))
			return null;
		for(int i = 0; i < column.length; i++)
			this.set(i, index, column[i]);
		return this;
	}
	public static MatrixForDic product(MatrixForDic arg0, MatrixForDic arg1) {
		if(arg0.getColumnDimension() != arg1.getRowDimension())
			return null;
		MatrixForDic result =
				new MatrixForDic(arg0.getRowDimension(), arg1.getColumnDimension());
		for(int i = 0; i < result.getRowDimension(); i++)
			for(int j = 0; j < result.getColumnDimension(); j++) {
				double element = 0;
				double[] row = arg0.getRowCopy(i),
						column = arg1.getColumnCopy(j);
				for(int k = 0; k < row.length; k++)
					element += row[k] * column[k];
				result.set(i, j, element);
			}
		return result;
	}
	public static double[] product(MatrixForDic arg0, double[] arg1) {
		if(arg0.getColumnDimension() != arg1.length)
			return null;
		double [] result = new double[arg0.getRowDimension()];
		for(int i = 0; i < result.length; i++) {
			double[] row = arg0.getRowCopy(i);
			double element = 0;
			for(int j = 0; j < row.length; j++)
				element += row[j] * arg1[j];
			result[i] = element;
		}
		return result;
	}
	public static MatrixForDic difference(MatrixForDic arg0, MatrixForDic arg1) {
		if((arg0.getRowDimension() != arg1.getRowDimension()) ||
				(arg0.getColumnDimension() != arg1.getColumnDimension()))
			return null;
		double[][] resultArray = arg0.getArrayCopy();
		for(int i = 0; i < arg1.getRowDimension(); i++)
			for(int j = 0; j < arg1.getColumnDimension(); j++)
				resultArray[i][j] -= arg1.get(i, j);
		return new MatrixForDic(resultArray);
	}
	public MatrixForDic copy() {
		return new MatrixForDic(this.getArrayCopy());
	}
	public MatrixForDic gInverse() {
		MatrixForDic t = new MatrixForDic(transpose().getArrayCopy());
		MatrixForDic inv = new MatrixForDic(product(t, this.copy()).inverse());
		return product(inv, t);
	}
	public MatrixForDic setMatrix(MatrixForDic matrix) {
		if((this.getRowDimension() != matrix.getRowDimension())
				|| (this.getColumnDimension() != matrix.getColumnDimension()))
			return null;
		for(int i = 0; i < this.getRowDimension(); i++)
			for(int j = 0; j < this.getColumnDimension(); j++)
				this.set(i, j, matrix.get(i, j));
		return this;
	}
	public MatrixForDic limitedWithColumn(int[] indexes) {
		MatrixForDic limited =
				new MatrixForDic(this.getRowDimension(), indexes.length);
		for(int i = 0; i < indexes.length; i++)
			if((indexes[i] < 0) || (indexes[i] <= getColumnDimension()))
				return null;
			else setColumn(i, this.getColumnCopy(indexes[i]));
		return limited;
	}
	public MatrixForDic limitedWidthRow(int[] indexes) {
		MatrixForDic limited =
				new MatrixForDic(indexes.length, this.getColumnDimension());
		for(int i = 0; i < indexes.length; i++)
			if((indexes[i] < 0) || (indexes[i] <= getRowDimension()))
				return null;
			else setRow(i, this.getRowCopy(indexes[i]));
		return limited;
	}
	public boolean isZero() {
		for(int i = 0; i < this.getRowDimension(); i++)
			for(int j = 0; j < this.getColumnDimension(); j++)
				if(this.get(i, j) != 0)
					return false;
		
		return true;
	}
}
