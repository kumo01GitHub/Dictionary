package monoimg;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import javax.imageio.ImageIO;

public class MonoImage {
	private BufferedImage image;
	private byte[][] monoArray;

	public MonoImage(String imgPath) {
		try {
			image = ImageIO.read(new File(imgPath));
		} catch (IOException e) {
			e.printStackTrace();
			image = null;
		}
		monoArray = toMonoArray(image);
	}
	public MonoImage(File imgFile) {
		try {
			image = ImageIO.read(imgFile);
		} catch (IOException e) {
			e.printStackTrace();
			image = null;
		}
		monoArray = toMonoArray(image);
	}
	public MonoImage(byte[][] array) {
		monoArray = array.clone();
		image = new BufferedImage(array.length, array[0].length, BufferedImage.TYPE_4BYTE_ABGR);
	}
	private byte[][] toMonoArray(BufferedImage img) {
		if(img == null) {
			return null;
		} else {
			byte[][] array = new byte[img.getWidth()][img.getHeight()];
			for(int i = 0; i < img.getWidth(); i++)
				for(int j = 0; j < img.getHeight(); j++) {
					Color color = new Color(img.getRGB(i, j));
					array[i][j] = ByteBuffer.allocate(4).putInt(
							(color.getRed() + color.getGreen() + color.getBlue()) / 3).get(3);
				}
			return array;
		}
	}
	public int getWidth() {
		if(monoArray == null)
			return -1;
		else
			return monoArray.length;
	}
	public int getHeight() {
		if(monoArray == null)
			return -1;
		else
			return monoArray[0].length;
	}
	public byte[][] getArrayCopy() {
		if(monoArray == null)
			return null;
		else
			return monoArray.clone();
	}
	public BufferedImage getMonoImage() {
		if(monoArray == null) {
			return null;
		} else {
			BufferedImage mono = new BufferedImage(
					image.getWidth(), image.getHeight(), image.getType());
			for(int i = 0; i < getWidth(); i++)
				for(int j = 0; j < getHeight(); j++) {
					byte[] pxByteArray = new byte[4];
					pxByteArray[0] = Byte.MAX_VALUE;
					for(int k = 1; k < 4; k++) pxByteArray[k] = monoArray[i][j];
					mono.setRGB(i, j, ByteBuffer.wrap(pxByteArray).getInt());
				}
			return mono;
		}
	}
	public File exportMonoImageAs(String name, String suffix) {
		BufferedImage mono = getMonoImage();
		if(mono == null) {
			return null;
		} else {
			File output = new File(name);
			try {
				ImageIO.write(mono, suffix, output);
			} catch (IOException e) {
				e.printStackTrace();
			}
			return output;
		}
	}
}
