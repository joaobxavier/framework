/*
 * Created on 23-okt-2003
 */
package nl.tudelft.bt.model.util;

import java.awt.Color;

/**
 * Implement colormaps
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class ColorMaps {
	private static final float GRAYVAL = 0.6f;
	private static final float BLACKVAL = 0.6f;

	/**
	 * Returns a color according to a float index in 0-1 range
	 * 
	 * @param f
	 */
	static public Color getJetColor(float f) {
		if (f > 0 && f <= 1) {
			float f2 = (f * 2.0f - 1.0f);
			float r = (f > 0.5f ? f2 : 0);
			float g = (f > 0.5f ? (1 - f2) : (1 + f2));
			float b = (f > 0.5f ? 0 : -f2);
			// normalize all by the highest value
			float n = Math.max(Math.max(r, g), b);
			return new Color(r / n, g / n, b / n);
		} else if (f == 0) {
			return Color.white;
		} else if (f < 0) {
			return Color.pink;
		} else if (f == Float.POSITIVE_INFINITY) {
			return Color.darkGray;
		} else if (f > 1) {
			// the saturation bove one
			return Color.black;
		}
		return Color.lightGray;
	}

	/**
	 * Returns a color according to a float index in 0-1 range
	 * 
	 * @param f
	 */
	static public Color getFullJetColor(float f) {
		if (f >= 0 && f <= 1) {
			float f2 = (f * 2.0f - 1.0f);
			float r = (f > 0.5f ? f2 : 0);
			float g = (f > 0.5f ? (1 - f2) : (1 + f2));
			float b = (f > 0.5f ? 0 : -f2);
			// normalize all by the highest value
			float n = Math.max(Math.max(r, g), b);
			return new Color(r / n, g / n, b / n);
		}
		return Color.lightGray;
	}

	/**
	 * Returns a color according to a float index in 0-1 range
	 * 
	 * @param f
	 */
	static public Color getBoneColor(float f) {
		float r = 1 - f * 0.8f;
		float g = 1 - ((float) Math.pow(f, 1.5)) * 0.8f;
		float b = g;
		// normalize all by the highest value
		return new Color(r, g, b);
	}

	/**
	 * Returns a color according to a float index in 0-1 range
	 * 
	 * @param f
	 * @param gamma
	 */
	static public Color getFullJetColor(float f, float gamma) {
		return getFullJetColor((float) Math.pow(f, gamma));
	}

	/**
	 * Returns a color according to a float index in 0-1 range from black (0) to
	 * red (1)
	 * 
	 * @param f
	 */
	static public Color getRedscaleColor(float f) {
		f = (f > 1 ? 1 : f);
		return new Color(f, 0, 0);
	}

	/**
	 * Returns a color according to a float index in 0-1 range
	 * 
	 * @param f
	 */
	static public Color getBluescaleColor(float f) {
		return new Color(0.5f * (1 - f), 0.5f * (1 - f), 1);
	}

	/**
	 * Return a color that is grayer than original color c by a factor f,
	 * keeping brightness and hue
	 * 
	 * @param c
	 * @param f
	 *            must be between 0 and 1
	 * @return color less saturated (grayer)
	 */
	static final public Color unsaturate(Color c, float f) {
		// make sure f is not negative nor greater than 1
		f = (f < 0 ? 0 : f);
		f = (f > 1 ? 1.0f : f);
		float r = (float) c.getRed() / 255.0f;
		float g = (float) c.getGreen() / 255.0f;
		float b = (float) c.getBlue() / 255.0f;
		return new Color((GRAYVAL - r) * f + r, (GRAYVAL - g) * f + g,
				(GRAYVAL - b) * f + b);
	}

	/**
	 * Reduce brightness of a color by factor f, keeping saturation and hue
	 * 
	 * @param c
	 * @param f
	 *            must be between 0 and 1
	 * @return darker color
	 */
	static final public Color darken(Color c, float f) {
		// make sure f is not negative nor greater than 1
		f = (f < 0 ? 0 : f);
		f = (f > 1 ? 1.0f : f);
		float r = (float) c.getRed() / 255.0f;
		float g = (float) c.getGreen() / 255.0f;
		float b = (float) c.getBlue() / 255.0f;
		return new Color(darkenValue(r, f), darkenValue(g, f),
				darkenValue(b, f));
	}

	static final private float darkenValue(float x, float f) {
		return ((BLACKVAL - 1.0f) * f + 1.0f) * x;
	}

	/**
	 * Increase brightness of a color by factor f, keeping saturation and hue
	 * 
	 * @param c
	 * @param f
	 *            must be between 0 and 1
	 * @return darker color
	 */
	static final public Color brighten(Color c, float f) {
		// make sure f is not negative nor greater than 1
		f = (f < 0 ? 0 : f);
		f = (f > 1 ? 1.0f : f);
		float r = (float) c.getRed() / 255.0f;
		float g = (float) c.getGreen() / 255.0f;
		float b = (float) c.getBlue() / 255.0f;
		return new Color(brightenValue(r, f), brightenValue(g, f),
				brightenValue(b, f));
	}

	static final public float brightenValue(float x, float f) {
		return x + (1 - x) * f;
	}

}
