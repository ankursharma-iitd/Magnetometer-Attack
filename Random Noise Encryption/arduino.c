#define RELAY1  7                        
void setup()

{    


Serial.begin(9600);
  pinMode(RELAY1, OUTPUT);       

}

  void loop()

{

   digitalWrite(RELAY1,0);           // Turns ON Relays 1
   Serial.println("Light ON");
   delay(random(1,2000));            // Random time delay

   digitalWrite(RELAY1,1);          // Turns Relay Off
   Serial.println("Light OFF");
   delay(random(1,2000));
   
}
