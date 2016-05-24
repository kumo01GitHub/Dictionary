package experiment;

import java.nio.ByteBuffer;

import dictionary.Dictionary;
import dictionary.MatrixForDic;
import monoimg.MonoImage;

public class Experiment {
	public static void main(String[] args) {
		MonoImage img = new MonoImage(args[0]);
		int dicSize = 32, sqSize = 8, sampleDim = 128;
		MatrixForDic samples = getRandomSamplesFromImg(img, sqSize, sampleDim);
		Dictionary dic = new Dictionary(samples, dicSize);
		System.out.println("making dictionary is done.");
		
		System.out.println("done.");
	}

	private static double[][] monoImgToArray(MonoImage img) {
		byte[][] bytes = img.getArrayCopy();
		double[][] array = new double[img.getWidth()][img.getHeight()];
		for(int i = 0; i < img.getWidth(); i++)
			for(int j = 0; j < img.getHeight(); j++)
				array[i][j] = (double)Byte.toUnsignedInt(bytes[i][j]);
		return array;
	}
	private static double[][] getSquareFromImg(MonoImage img, int x, int y, int size) {
		if((img.getWidth() < x + size) || (img.getHeight() < y + size)) {
			return null;
		} else {
			double[][] square = new double[size][size],
					imgArray = monoImgToArray(img);
			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
					square[i][j] = imgArray[x + i][y + j];
			return square;
		}
	}
	private static double[] array2Dto1D(double[][] array) {
		if(array == null) {
			return null;
		} else {
			double[] array1 =
					new double[array.length * array[0].length];
			for(int i = 0; i < array.length; i++)
				for(int j = 0; j < array[0].length; j++)
					array1[array[0].length * i + j] = array[i][j];
			return array1;
		}
	}
	private static double[][] array1Dto2D(double[] array, int width, int height) {
		if(array == null) {
			return null;
		} else if(array.length != width * height) {
			return null;
		} else {
			double[][] array2 = new double[width][height];
			for(int i = 0; i < width; i++)
				for(int j = 0; j < height; j++)
					array2[i][j] = array[i * height + j];
			return array2;
		}
	}
	private static MatrixForDic getAllSamplesFromImg(MonoImage img, int size) {
		if(img == null) {
			return null;
		} else if((img.getWidth() < size) || (img.getHeight() < size)) {
			return null;
		} else {
			MatrixForDic samples =
					new MatrixForDic(size*size, (img.getWidth()-size+1)*(img.getHeight()-size+1));
			for(int x = 0; x < img.getWidth()-size+1; x++)
				for(int y = 0; y < img.getHeight()-size+1; y++)
					samples.setColumn((img.getHeight()-size+1)*x + y, array2Dto1D(getSquareFromImg(img, x, y, size)));
			return samples;
		}
	}
	private static MatrixForDic getRandomSamplesFromImg(MonoImage img, int size, int dim) {
		if(img == null) {
			return null;
		} else if((img.getWidth() < size) || (img.getHeight() < size)) {
			return null;
		} else {
			MatrixForDic samples = new MatrixForDic(size*size, dim);
			for(int i = 0; i < dim; i++)
				samples.setColumn(i, array2Dto1D(getSquareFromImg(img,
						(int)(Math.random()*(img.getWidth() - size + 1)),
						(int)(Math.random()*(img.getHeight() - size + 1)), size)));
			return samples;
		}
	}
	private static double maxAbsElement(MatrixForDic matrix) {
		if(matrix == null) {
			return -1;
		} else {
			double maxAbs = 0;
			for(int i = 0; i < matrix.getRowDimension(); i++)
				for(int j = 0; j < matrix.getColumnDimension(); j++) {
					double abs = Math.abs(matrix.get(i, j));
					if(maxAbs < abs) maxAbs = abs;
				}
			return maxAbs;
		}
	}
	private static MonoImage arrayToMonoImg(double[][] array) {
		if(array == null) {
			return null;
		} else if(array[0] == null) {
			return null;
		} else {
			byte[][] monoArray = new byte[array.length][array[0].length];
			for(int i = 0; i < monoArray.length; i++)
				for(int j = 0; j < monoArray[0].length; j++)
					monoArray[i][j] = ByteBuffer.allocate(4).putInt((int)array[i][j]).get(3);
			return new MonoImage(monoArray);
		}
	}
}
