/* parameters:
 * 	omp.maxWT
 * 	omp.error
 *	makeDictionary.error
 *	makeDictionary.repeatMax
 *	makeSparseMatrix.error
 *	makeSparseMatrix.repeatMax
*/

package dictionary;

import java.util.ArrayList;
import java.util.List;
import Jama.SingularValueDecomposition;

public class Dictionary extends MatrixForDic {
	private MatrixForDic originalSamples, originalSpMatrix;

	public Dictionary(MatrixForDic samples, int size) {
		super(samples.getRowDimension(), size);
		makeDictionary(samples);
	}

	public int getDicSize() {
		return this.getColumnDimension();
	}
	private MatrixForDic omp(MatrixForDic spMatrix) {
		int maxWT = 5;
		double error = 0.1;
		MatrixForDic updatedSpMatrix =
				new MatrixForDic(spMatrix.getRowDimension(), spMatrix.getColumnDimension());
		for(int sampleIDX = 0; sampleIDX < originalSamples.getColumnDimension(); sampleIDX++) {
			double[] sample = originalSamples.getColumnCopy(sampleIDX);
			List<Integer> spIDXes = new ArrayList<Integer>();
			double[] updatedCol = new double[getDicSize()];
			for(int i = 0; i < updatedCol.length; i++) updatedCol[i] = 0;
			for(int wt = 0; wt < maxWT; wt++) {
				double[] errorArray = sample.clone();
				for(int i = 0; i < errorArray.length; i++)
					errorArray[i] -= MatrixForDic.product(this, updatedCol)[i];
				double minD = Double.MAX_VALUE;
				for(int i = 0; i < getDicSize(); i++) {
					if(spIDXes.contains(i)) continue;
					double errNorm = 0, dicColNorm = 0, pDicErr = 0;
					for(int l = 0; l < errorArray.length; l++) {
						errNorm += Math.pow(errorArray[l], 2);
						dicColNorm += Math.pow(this.get(l, i), 2);
						pDicErr += (this.get(l, i) * errorArray[l]);
					}
					double d = (errNorm * dicColNorm) - (pDicErr * pDicErr);
					if(d < minD) {
						if(wt < spIDXes.size()) spIDXes.remove(spIDXes.size() - 1);
						spIDXes.add(i);
						minD = d;
					}
				}
				MatrixForDic limDic = new MatrixForDic(this.getRowDimension(), spIDXes.size());
				for(int i = 0; i < limDic.getColumnDimension(); i++)
					limDic.setColumn(i, this.getColumnCopy(spIDXes.get(i)));
				MatrixForDic limDicInv = limDic.gInverse();
				double[] tmpUpdatedCol = MatrixForDic.product(limDicInv, sample);
				for(int i = 0; i < updatedCol.length; i++) updatedCol[i] = 0;
				for(int i = 0; i < spIDXes.size(); i++)
					updatedCol[spIDXes.get(i)] = tmpUpdatedCol[i];
				double tmpError = 0;
				for(int i = 0; i < sample.length; i++)
					tmpError += Math.pow(sample[i] - MatrixForDic.product(this, updatedCol)[i], 2);
				if(tmpError < Math.pow(error, 2)) break;
			}
			updatedSpMatrix.setColumn(sampleIDX, updatedCol);
		}
		
		return updatedSpMatrix;
	}
	private MatrixForDic omp(MatrixForDic spMatrix, MatrixForDic samples) {
		int maxWT = 5;
		double error = 0.1;
		MatrixForDic updatedSpMatrix =
				new MatrixForDic(spMatrix.getRowDimension(), spMatrix.getColumnDimension());
		for(int sampleIDX = 0; sampleIDX < samples.getColumnDimension(); sampleIDX++) {
			double[] sample = samples.getColumnCopy(sampleIDX);
			List<Integer> spIDXes = new ArrayList<Integer>();
			double[] updatedCol = new double[getDicSize()];
			for(int i = 0; i < updatedCol.length; i++) updatedCol[i] = 0;
			for(int wt = 0; wt < maxWT; wt++) {
				double[] errorArray = sample.clone();
				for(int i = 0; i < errorArray.length; i++)
					errorArray[i] -= MatrixForDic.product(this, updatedCol)[i];
				double minD = Double.MAX_VALUE;
				for(int i = 0; i < getDicSize(); i++) {
					if(spIDXes.contains(i)) continue;
					double errNorm = 0, dicColNorm = 0, pDicErr = 0;
					for(int l = 0; l < errorArray.length; l++) {
						errNorm += Math.pow(errorArray[l], 2);
						dicColNorm += Math.pow(this.get(l, i), 2);
						pDicErr += (this.get(l, i) * errorArray[l]);
					}
					double d = (errNorm * dicColNorm) - (pDicErr * pDicErr);
					if(d < minD) {
						if(wt < spIDXes.size()) spIDXes.remove(spIDXes.size() - 1);
						spIDXes.add(i);
						minD = d;
					}
				}
				MatrixForDic limDic = new MatrixForDic(this.getRowDimension(), spIDXes.size());
				for(int i = 0; i < limDic.getColumnDimension(); i++)
					limDic.setColumn(i, this.getColumnCopy(spIDXes.get(i)));
				MatrixForDic limDicInv = limDic.gInverse();
				double[] tmpUpdatedCol = MatrixForDic.product(limDicInv, sample);
				for(int i = 0; i < updatedCol.length; i++) updatedCol[i] = 0;
				for(int i = 0; i < spIDXes.size(); i++)
					updatedCol[spIDXes.get(i)] = tmpUpdatedCol[i];
				double tmpError = 0;
				for(int i = 0; i < sample.length; i++)
					tmpError += Math.pow(sample[i] - MatrixForDic.product(this, updatedCol)[i], 2);
				if(tmpError < Math.pow(error, 2)) break;
			}
			updatedSpMatrix.setColumn(sampleIDX, updatedCol);
		}
		
		return updatedSpMatrix;
	}
	private MatrixForDic ksvd(MatrixForDic spMatrix) {
		MatrixForDic updatedDic = this.copy(),
				updatedSp = new MatrixForDic(spMatrix.getRowDimension(), spMatrix.getColumnDimension());
		for(int k = 0; k < this.getDicSize(); k++) {
			List<Integer> spIDXes = new ArrayList<Integer>();
			for(int i = 0; i < spMatrix.getColumnDimension(); i++)
				if(spMatrix.get(k, i) != 0) spIDXes.add(i);
			if(spIDXes.size() == 0) continue;
			MatrixForDic lostDic =
					new MatrixForDic(this.getRowDimension(), this.getDicSize() - 1),
					lostSp = new MatrixForDic(spMatrix.getRowDimension() - 1, spMatrix.getColumnDimension());
			for(int i = 0; i < this.getDicSize() - 1; i++)
				if(i < k) {
					lostDic.setColumn(i, this.getColumnCopy(i));
					lostSp.setRow(i, spMatrix.getRowCopy(i));
				} else {
					lostDic.setColumn(i, this.getColumnCopy(i + 1));
					lostSp.setRow(i, spMatrix.getRowCopy(i + 1));
				}
			MatrixForDic lostDicSp = MatrixForDic.product(lostDic, lostSp),
					limLostDicSp = new MatrixForDic(lostDicSp.getRowDimension(), spIDXes.size()),
					limSamples = new MatrixForDic(originalSamples.getRowDimension(), spIDXes.size());
			for(int i = 0; i < spIDXes.size(); i++) {
				limLostDicSp.setColumn(i, lostDicSp.getColumnCopy(spIDXes.get(i)));
				limSamples.setColumn(i, originalSamples.getColumnCopy(spIDXes.get(i)));
			}
			MatrixForDic limErr = MatrixForDic.difference(limSamples, limLostDicSp);
			MatrixForDic u, v;
			double sv;
			if(limErr.getRowDimension() < limErr.getColumnDimension()) {
				limErr = new MatrixForDic(limErr.transpose());
				SingularValueDecomposition svd = limErr.svd();
				u = new MatrixForDic((svd.getV()).transpose());
				v = new MatrixForDic(svd.getU());
				sv = (svd.getS()).get(0, 0);
			} else {
				SingularValueDecomposition svd = limErr.svd();
				u = new MatrixForDic(svd.getU());
				v = new MatrixForDic(svd.getV());
				sv = (svd.getS()).get(0, 0);
			}
			updatedDic.setColumn(k, u.getColumnCopy(0));
			double[] updatedSpRow = spMatrix.getRowCopy(k);
			double[] vCol = v.getColumnCopy(0);
			for(int i = 0; i < spIDXes.size(); i++)
				updatedSpRow[spIDXes.get(i)] = sv * vCol[i];
			updatedSp.setRow(k, updatedSpRow);
		}
		this.setMatrix(updatedDic);
		return updatedSp;
	}
	private MatrixForDic ksvd(MatrixForDic spMatrix, MatrixForDic samples) {
		MatrixForDic updatedSp = new MatrixForDic(spMatrix.getRowDimension(), spMatrix.getColumnDimension());
		for(int k = 0; k < this.getDicSize(); k++) {
			List<Integer> spIDXes = new ArrayList<Integer>();
			for(int i = 0; i < spMatrix.getColumnDimension(); i++)
				if(spMatrix.get(k, i) != 0) spIDXes.add(i);
			if(spIDXes.size() == 0) continue;
			MatrixForDic lostDic =
					new MatrixForDic(this.getRowDimension(), this.getDicSize() - 1),
					lostSp = new MatrixForDic(spMatrix.getRowDimension() - 1, spMatrix.getColumnDimension());
			for(int i = 0; i < this.getDicSize() - 1; i++)
				if(i < k) {
					lostDic.setColumn(i, this.getColumnCopy(i));
					lostSp.setRow(i, spMatrix.getRowCopy(i));
				} else {
					lostDic.setColumn(i, this.getColumnCopy(i + 1));
					lostSp.setRow(i, spMatrix.getRowCopy(i + 1));
				}
			MatrixForDic lostDicSp = MatrixForDic.product(lostDic, lostSp),
					limLostDicSp = new MatrixForDic(lostDicSp.getRowDimension(), spIDXes.size()),
					limSamples = new MatrixForDic(samples.getRowDimension(), spIDXes.size());
			for(int i = 0; i < spIDXes.size(); i++) {
				limLostDicSp.setColumn(i, lostDicSp.getColumnCopy(spIDXes.get(i)));
				limSamples.setColumn(i, samples.getColumnCopy(spIDXes.get(i)));
			}
			MatrixForDic limErr = MatrixForDic.difference(limSamples, limLostDicSp);
			MatrixForDic v;
			double sv;
			if(limErr.getRowDimension() < limErr.getColumnDimension()) {
				limErr = new MatrixForDic(limErr.transpose());
				SingularValueDecomposition svd = limErr.svd();
				v = new MatrixForDic(svd.getU());
				sv = (svd.getS()).get(0, 0);
			} else {
				SingularValueDecomposition svd = limErr.svd();
				v = new MatrixForDic(svd.getV());
				sv = (svd.getS()).get(0, 0);
			}
			double[] updatedSpRow = spMatrix.getRowCopy(k);
			double[] vCol = v.getColumnCopy(0);
			for(int i = 0; i < spIDXes.size(); i++)
				updatedSpRow[spIDXes.get(i)] = sv * vCol[i];
			updatedSp.setRow(k, updatedSpRow);
		}
		return updatedSp;
	}
	private Dictionary makeDictionary(MatrixForDic samples) {
		double error = 0.1;
		int repeatMax = 128;
		originalSamples = new MatrixForDic(samples.getArrayCopy());
		originalSpMatrix =
				new MatrixForDic(getDicSize(), samples.getColumnDimension());
		randomInit();
		for(int i = 0; i < repeatMax; i++) {
			originalSpMatrix = omp(originalSpMatrix);
			originalSpMatrix = ksvd(originalSpMatrix);
			if((MatrixForDic.difference(originalSamples, MatrixForDic.product(this, originalSpMatrix))).norm2() < error)
				break;
		}
		return this;
	}
	private void randomInit() {
		for(int i = 0; i < this.getColumnDimension(); i++) {
			double[] col = new double[this.getRowDimension()];
			for(int j = 0; j < col.length; j++)
				col[j] = Math.random();
			
			double norm = 0;
			for(int j = 0; j < col.length; j++)
				norm += Math.pow(col[j], 2);
			norm = Math.pow(norm, 0.5);
			
			for(int j = 0; j < col.length; j++)
				col[j] /= norm;
			this.setColumn(i, col);
		}
	}
	public MatrixForDic getSparseMatrixCopy() {
		return originalSpMatrix.copy();
	}
	public MatrixForDic makeSparseMatrix(MatrixForDic samples) {
		double error = 0.1;
		int repeatMax = 128;
		MatrixForDic spMatrix = new MatrixForDic(getDicSize(), samples.getColumnDimension());
		for(int i = 0; i < repeatMax; i++) {
			spMatrix = omp(spMatrix, samples);
			spMatrix = ksvd(spMatrix, samples);
			if((MatrixForDic.difference(samples, MatrixForDic.product(this, spMatrix))).norm2() < error)
				break;
		}
		return spMatrix;
	}
}
