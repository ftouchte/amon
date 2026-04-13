package io.github.ftouchte.fitting;

public class FitParameter {
    double value;
    double minValue;
    double maxValue;

    public FitParameter(double _value, double _min, double _max) {
        value = _value;
        minValue = _min;
        maxValue = _max;
    }

    public double getValue() {
        return value;
    }

    public double getMinValue() {
        return minValue;
    }

    public double getMaxValue() {
        return maxValue;
    }

    public void setParameter(double _value, double _min, double _max) {
        value = _value;
        minValue = _min;
        maxValue = _max;
    }

    String name = "";

    public void setName(String _name) {
        name = _name;
    }

    public String getName() {
        return name;
    }

    public void print() {
        System.out.printf("  (%s) value : %f , min : %f , max : %f\n", name.equals("") ? "unnammed" : name, value, minValue, maxValue);
    }
}
