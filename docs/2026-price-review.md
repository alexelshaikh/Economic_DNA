# 2026 price baseline and decline assumptions

Reviewed 2026-09-27. The website baseline has been updated; the paper's
`storage_service.py` reference implementation has not yet been rebased.

## Applied price updates

- Cloud and tape price base years and the default scenario start year: 2026.
- AWS and Azure annual decline defaults remain 10%; these are scenario assumptions,
  not observed 2025-to-2026 price reductions.
- AWS: us-east-1, Deep Archive, bulk restoration. Azure: East US 2,
  Archive LRS, flat namespace, standard-priority retrieval.
- Cloud capacity rates use USD per decimal MB, converted from provider binary GB
  using `USD_per_GiB * 1_000_000 / 2**30`. Request charges are unchanged.
- Tape media: HPE Q2079A, listed by BlueAlly on HPE Store at USD 104.14,
  divided by 18 TB native capacity: USD 5.785555555555556/TB.
  This is a dated reseller quote, excluding shipping, taxes and additional copies.
- Tape hardware and energy prices remain model assumptions, not refreshed quotes.
- DNA prices and fits are unchanged.
- Currency metadata identifies USD and a 2026 price baseline, without claiming
  constant purchasing power or a PPP adjustment. No inflation conversion was made.

Price sources:

- [AWS regional price catalog](https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonS3GlacierDeepArchive/current/us-east-1/index.json)
- [AWS billing units](https://aws.amazon.com/s3/pricing/)
- [Azure retail prices](https://prices.azure.com/api/retail/prices?$filter=armRegionName%20eq%20%27eastus2%27%20and%20skuName%20eq%20%27Archive%20LRS%27%20and%20productName%20eq%20%27Azure%20Data%20Lake%20Storage%20Gen2%20Flat%20Namespace%27)
- [Azure billing units](https://azure.microsoft.com/en-us/pricing/details/storage/blobs/)
- [HPE media listing](https://buy.hpe.com/us/en/storage/storage-media/tape-media/hpe-storage-lto-ultrium-cartridges/p/34648)
- [HPE native capacity specification](https://www.hpe.com/us/en/collaterals/collateral.c04154496.html)

## Tape decline recommendation (not yet applied)

Recommend media 10%, hardware 0%, energy 0% per year as explicit modeling
assumptions, not three empirically estimated rates. Existing website tape
declines remain 20%, 0%, 15% pending agreement.

Lauhoff and Shahidi, IBM, *Data Storage Trends: NAND, HDD and Tape Storage*,
Library of Congress, March 9, 2026, slide 14 reports a 14% annual LTO media
cost/bit decline over 2012-2025. Slides 9 and 16 qualify the series: lowest
available cartridge prices, compressed capacity, and no drives or libraries.
Thus 10% is a slower-improvement scenario, not the reported estimate. A native
capacity trend would require accounting for compression assumptions over time.
[Presentation](https://digitalpreservation.gov/meetings/DSA2026/0106_LoC2026_IBM_Landscape_GL.pdf)

Hardware: no comparable verified time series of library/drive cost per usable TB
was established in this review. Holding its allocated cost fixed is a transparent
modeling choice, not evidence that equipment prices cannot change.

Energy: EIA's September 2026 Short-Term Energy Outlook, Table 7a, forecasts US
commercial electricity at 13.99 cents/kWh in 2026 and 14.15 in 2027, versus
13.41 in 2025. These are near-term nominal tariffs, not tape energy per TB.
Energy cost/TB/year also depends on kWh/TB/year and utilization; no source here
establishes a 15% annual decline in that product. A 0% net change is a simplifying
baseline, not an EIA forecast for tape.
[EIA report](https://www.eia.gov/outlooks/steo/pdf/steo_full.pdf)

For the paper, test media decline rates of 0%, 10%, 14%, and 20%, clearly labeled
as chosen sensitivity cases. These are not confidence bounds. Century-scale
exponential decline should be framed as an extrapolation, not a forecast;
include plateau or finite-duration improvement cases before making robust
long-horizon claims. Higher conventional-storage decline generally disadvantages
DNA, so 'conservative' should always identify which technology it favors.
