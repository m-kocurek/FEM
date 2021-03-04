

public class Node {
        private double x;
        private double y;
        private boolean status;
        private final int uid;
        private double temp;

        public Node(double x, double y, final int uid , boolean status, double initialTemp) {
            this.x = x;
            this.y = y;
            this.uid = uid;
            this.status = status;
            this.temp = initialTemp;
        }

        public Node(Point point){
            this.x = point.getX();
            this.y = point.getY();
            uid = -1;
            status = false;
        }

        public double getTemp() {
            return temp;
        }

        public void setTemp(double temp) {
            this.temp = temp;
        }

        public boolean isStatus() {
            return status;
        }

        public double getX() {
            return x;
        }
        public double getY() {
            return y;
        }
}
