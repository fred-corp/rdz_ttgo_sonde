
/* Horus decoder functions */
#include "Horus.h"
#include "SX1278FSK.h"
#include "rsc.h"
#include "Sonde.h"

#define Horus_DEBUG 1

#if Horus_DEBUG
#define Horus_DBG(x) x
#else
#define Horus_DBG(x)
#endif

#define HorusMAXLEN (128)
static byte data[800];
static int dpos = 0;


// whole 51 row frame as C structure
// taken from https://github.com/einergehtnochrein/ra-firmware
struct subframeBuffer {
    uint64_t valid;   // bitmask for subframe valid; lsb=frame 0, etc.
    union {
        byte rawData[51*16];
        struct __attribute__((__packed__)) {
            uint16_t crc16;                     /* CRC16 CCITT Checksum over range 0x002...0x31F */
            uint16_t frequency;                 /* 0x002: TX is on 400 MHz + (frequency / 64) * 10 kHz */
            uint8_t startupTxPower;             /* 0x004: TX power level at startup (1...7) */
            uint8_t reserved005;
            uint8_t reserved006;
            uint16_t reserved007;               /* 0x007:  ?? (some bitfield) [0],[1],[2],[3]. Init value = 0xE */
            uint16_t reserved009;               /* 0x009: ? */
            uint8_t reserved00B;
            uint8_t reserved00C;
            uint8_t serial[8];                  /* 0x00D: Sonde ID, 8 char, not terminated */
            uint16_t firmwareVersion;           /* 0x015: 10000*major + 100*minor + patch*/
            uint16_t reserved017;
            uint16_t minHeight4Flight;          /* 0x019: Height (meter above ground) where flight mode begins */
            uint8_t lowBatteryThreshold100mV;   /* 0x01B: (Default=18) Shutdown if battery voltage below this
                                                          threshold for some time (10s ?)
                                                */
            uint8_t nfcDetectorThreshold;       /* 0x01C: NFC detector threshold [25mV] (Default: 0x05 = 125mV) */
            uint8_t reserved01D;                /* 0x01D: ?? (Init value = 0xB4) */
            uint8_t reserved01E;                /* 0x01E: ?? (Init value = 0x3C) */
            uint16_t reserved01F;
            int8_t refTemperatureThreshold;     /* 0x021: Reference temperature threshold [°C] */
            uint8_t reserved022;
            uint16_t reserved023;
            uint16_t reserved025;
            int16_t flightKillFrames;           /* 0x027: Number of frames in flight until kill (-1 = disabled) */
            uint16_t reserved029;               /* 0x029: ? (Init value = 0) */
            uint8_t burstKill;                  /* 0x02B: Burst kill (0=disabled, 1=enabled) */
            uint8_t reserved02C;
            uint8_t reserved02D;
            uint16_t reserved02E;
            uint16_t reserved030;
            uint8_t reserved032;
            uint16_t reserved033;
            uint16_t reserved035;
            uint16_t reserved037;
            uint16_t reserved039;               /* 0x039: */
            uint8_t reserved03B;                /* 0x03B: */
            uint8_t reserved03C;                /* 0x03C: */
            float refResistorLow;               /* 0x03D: Reference resistor low (750 Ohms) */
            float refResistorHigh;              /* 0x041: Reference resistor high (1100 Ohms) */
            float refCapLow;                    /* 0x045: Reference capacitance low (0) */
            float refCapHigh;                   /* 0x049: Reference capacitance high (47 pF) */
            float taylorT[3];                   /* 0x04D: Tayor coefficients for main temperature calculation */
            float calT;                         /* 0x059: Calibration factor for main sensor */
            float polyT[6];                     /* 0x05D: */
            float calibU[2];                    /* 0x075: Calibration coefficients for humidity sensor */

            float matrixU[7][6];                /* 0x07D: Matrix for humidity sensor RH calculation */
            float taylorTU[3];                  /* 0x125: Coefficients for U sensor temperature calculation */
            float calTU;                        /* 0x131: Calibration factor for U temperature sensor */
            float polyTrh[6];                   /* 0x135:  */

            uint8_t reserved14D;                /* 0x14D: */
            uint32_t reserved14E;               /* 0x14E: */

            float f152;
            uint8_t u156;
            float f157;                         /* 0x157: ?? (Initialized by same value as calibU[0]) */
            uint8_t reserved15B;                /* 0x15B: */
            uint32_t reserved15C;               /* 0x15C: */
            float f160[35];
            uint8_t startIWDG;                  /* 0x1EC: If ==1 or ==2: Watchdog IWDG will not be started */
            uint8_t parameterSetupDone;         /* 0x1ED: Set (!=0) if parameter setup was done */
            uint8_t enableTestMode;             /* 0x1EE: Test mode (service menu) (0=disabled, 1=enabled) */
            uint8_t enableTX;                   /* 0x1EF: 0=TX disabled, 1=TX enabled (maybe this is autostart?) */
            float f1F0[8];
            float pressureLaunchSite[2];        /* 0x210: Pressure [hPa] at launch site */
            struct __attribute__((__packed__)){
                char variant[10];               /* 0x218: Sonde variant (e.g. "Horus-SG") */
                uint8_t mainboard[10];          /* 0x222: Name of mainboard (e.g. "RSM412") */
            } names;
            struct __attribute__((__packed__)){
                uint8_t mainboard[9];           /* 0x22C: Serial number of mainboard (e.g. "L1123553") */
                uint8_t text235[12];            /* 0x235: "0000000000" */
                uint16_t reserved241;           /* 0x241: */
                uint8_t pressureSensor[8];      /* 0x243: Serial number of pressure sensor (e.g. "N1310487") */
                uint16_t reserved24B;           /* 0x24B: */
            } serials;
            uint16_t reserved24D;               /* 0x24D: */
            uint16_t reserved24F;               /* 0x24F: */
            uint16_t reserved251;               /* 0x251: (Init value = 0x21A = 538) */
            uint8_t xdataUartBaud;              /* 0x253: 1=9k6, 2=19k2, 3=38k4, 4=57k6, 5=115k2 */
            uint8_t reserved254;
            float cpuTempSensorVoltageAt25deg;  /* 0x255: CPU temperature sensor voltage at 25°C */
            uint8_t reserved259;
            uint8_t reserved25A[0x25E -0x25A];
            float matrixP[18];                  /* 0x25E: Coefficients for pressure sensor polynomial */
            float vectorBp[3];                  /* 0x2A6: */
            uint8_t reserved2B2[8];             /* 0x2B2: */
            float matrixBt[12];                 /* 0x2BA: */
            uint8_t reserved2EA[0x2FA-0x2EA];
            uint16_t halfword2FA[9];
            float reserved30C;
            float reserved310;                  /* 0x310: */
            uint8_t reserved314;                /* 0x314: */
            uint8_t reserved315;                /* 0x315: */
            int16_t burstKillFrames;            /* 0x316: Number of active frames after burst kill */
            uint8_t reserved318[0x320-0x318];

            /* This is fragment 50. It only uses 14 valid bytes! */
            int16_t killCountdown;              /* 0x320: Counts frames remaining until kill (-1 = inactive) */
            uint8_t reserved322[6];
            int8_t intTemperatureCpu;           /* 0x328: Temperature [°C] of CPU */
            int8_t intTemperatureRadio;         /* 0x329: Temperature [°C] of radio chip */
            int8_t reserved32A;                 /* 0x32A: */
            uint8_t reserved32B;                /* 0x32B: */
            uint8_t reserved32C;                /* 0x32C: ? (the sum of two slow 8-bit counters) */
            uint8_t reserved32D;                /* 0x32D: ? (the sum of two slow 8-bit counters) */
        } value;
    };
};
// moved global variable "calibration" to sondeInfo->extra

static uint16_t CRCTAB[256];

#define X2C_DIVR(a, b) ((b) != 0.0f ? (a)/(b) : (a))
#define X2C_DIVL(a, b) ((a)/(b))
static uint32_t X2C_LSH(uint32_t a, int32_t length, int32_t n)
{
	uint32_t m;

	m = 0;
	m = (length == 32) ? 0xFFFFFFFFl : (1 << length) - 1;
	if (n > 0) {
		if (n >= (int32_t)length)
			return 0;
		return (a << n) & m;
	}

	if (n <= (int32_t)-length)
		return 0;
	return (a >> -n) & m;
}


static void Gencrctab(void)
{
   uint16_t j;
   uint16_t i;
   uint16_t crc;
   for (i = 0U; i<=255U; i++) {
      crc = (uint16_t)(i*256U);
      for (j = 0U; j<=7U; j++) {
         if ((0x8000U & crc)) crc = X2C_LSH(crc,16,1)^0x1021U;
         else crc = X2C_LSH(crc,16,1);
      } /* end for */
      CRCTAB[i] = X2C_LSH(crc,16,-8)|X2C_LSH(crc,16,8);
   } /* end for */
} /* end Gencrctab() */

decoderSetupCfg horusSetupCfg = {
	.bitrate = 100,
	.rx_cfg = 0x1E, // Enable auto-AFC, auto-AGC, RX Trigger by preamble
	.sync_cfg = 0x57, // Set autostart_RX to 01, preamble 0, SYNC detect==on, syncsize=3 (==4 byte
	.sync_len = 8,
	.sync_data = (const uint8_t *)"\x08\x6D\x53\x88\x44\x69\x48\x1F",
	.preamble_cfg = 0xA8,
};


int Horus::setup(float frequency, int /*type*/) 
{
	if(!initialized) {
		Gencrctab();
		initrsc();
		initialized = true;
	}

	if(sx1278.ON()!=0) {
		Horus_DBG(Serial.println("Setting SX1278 power on FAILED"));
		return 1;
	}
	if(DecoderBase::setup(horusSetupCfg, sonde.config.horus.agcbw, sonde.config.horus.rxbw)!=0 ) {
		return 1;
	}

	// Packet config 1: fixed len, no mancecer, no crc, no address filter
	// Packet config 2: packet mode, no home ctrl, no beackn, msb(packetlen)=0)
	if(sx1278.setPacketConfig(0x08, 0x40)!=0) {
		Horus_DBG(Serial.println("Setting Packet config FAILED"));
		return 1;
	}
	int retval = sx1278.setFrequency(frequency);
	dpos = 0;

        sx1278.clearIRQFlags();

	return retval;
}

uint32_t Horus::bits2val(const uint8_t *bits, int len) {
	uint32_t val = 0;
	for (int j = 0; j < len; j++) {
		val |= (bits[j] << (len-1-j));
	}
	return val;
}

Horus::Horus() {
}


void Horus::printRaw(uint8_t *data, int len)
{
	char buf[3];
	int i;
	for(i=0; i<len; i++) {
		snprintf(buf, 3, "%02X", data[i]);
		Serial.print(buf);
	}
	Serial.println();
}

void Horus::bitsToBytes(uint8_t *bits, uint8_t *bytes, int len)
{
	int i;
	for(i=0; i<len*4; i++) {
	       	bytes[i/8] = (bytes[i/8]<<1) | (bits[i]?1:0);
	}
	bytes[(i-1)/8] &= 0x0F;
}


int Horus::receive() {
	sx1278.setPayloadLength(HorusMAXLEN-8); 
	int e = sx1278.receivePacketTimeout(5000, data+8);
#if 0
	if(e) { /*Serial.println("TIMEOUT");*/ return RX_TIMEOUT; } 

        for(int i=0; i<HorusMAXLEN; i++) { data[i] = reverse(data[i]); }
        for(int i=0; i<HorusMAXLEN; i++) { data[i] = data[i] ^ scramble[i&0x3F]; }
        Serial.print("printRaw : ");
        printRaw(data, HorusMAXLEN);
        return 0;
#else
  for(int i=0; i<HorusMAXLEN; i++) { 
    Serial.print(data[i], HEX);
    i % 8 > 0 ? Serial.print(" ") : Serial.println();
  }
  Serial.println();
	// FAKE testing data
	SondeInfo *si = sonde.si();
	si->d.lat = 50.78;
	si->d.lon = 4.77;
	si->d.alt = 250;
  si->d.hs = 2;
	si->d.vs = 3.4;
	si->d.validPos = 0x7f;
	si->d.validID = 1;
	strcpy(si->d.id, "Test!");
	return 0;
#endif
}

int Horus::waitRXcomplete() {
	// Currently not used. can be used for additinoal post-processing
	// (required for RS92 to avoid FIFO overrun in rx task)
	return 0;
}

// copy variant string to buf (max buflen chars; buflen should be 11
// return 0 if subtype is available, -1 if not
int Horus::getSubtype(char *buf, int buflen, SondeInfo *si) {
	struct subframeBuffer *sf = (struct subframeBuffer *)si->extra;
	if(!sf) return -1;
	if( ( (sf->valid>>0x21) &3) != 3 ) return -1;   // or 1 instead of 3 for the first 8 chars only, as in autorx?
	if(buflen>11) buflen=11;			    // then buflen should be capped at 9 (8+trailing \0)
	strncpy(buf, sf->value.names.variant, buflen);
	buf[buflen-1]=0;
	if(*buf==0) return -1;
	Serial.printf("subframe valid: %x%08x; subtype=%s\n", (uint32_t)(sf->valid>>32), (uint32_t)sf->valid, buf);
	return 0;
}

Horus horus = Horus();
