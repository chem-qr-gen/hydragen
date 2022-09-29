from pymongo import MongoClient

from mongo_address import mongo_address # secret address to mongodb database, change or remove for local testing

client = MongoClient(mongo_address)
db = client.chemquest_db
ms_qns = db.ms_qns

questions = [
    {
        "qid": 1,
        "text1": "Find the gaseous compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/XtX944c.gif",
        "text2": "",
        "hints": ["The compound does not contain carbon or hydrogen.", "The compound is a common contributor to acid rain."],
        "answers": ["sulfur dioxide", "SO2"]
    },
    {
        "qid": 2,
        "text1": "Find the gaseous compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/3okdkqT.gif",
        "text2": "",
        "hints": ["How many hydrogens are there in this compound?", "The compound is used to preserve tissues and organs."],
        "answers": ["formaldehyde", "methanal", "CH2O"]
    },
    {
        "qid": 3,
        "text1": "Find the pure liquid compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/M4o462D.gif",
        "text2": "",
        "hints": ["The compound is symmetrical about its centre.", "What do the peaks at m/z=15 and m/z=43 commonly represent?"],
        "answers": ["2,3-butadione", "buta-2,3-dione", "CH3COCOCH3"]
    },
    {
        "qid": 4,
        "text1": "Find the pure liquid compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/SMrzZgJ.gif",
        "text2": "",
        "hints": ["m/z=30 and m/z=46 are signs of what functional group?", "The compound is used as a fuel additive in motor sports."],
        "answers": ["nitromethane", "CH3NO2"]
    },
    {
        "qid": 5,
        "text1": "Find the pure liquid compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/3G9TqYP.gif",
        "text2": "",
        "hints": ["The compound is aromatic.", "What does the m/z of the molecular ion tell you about what atoms are present?"],
        "answers": ["pyridine", "C5H5N"]
    },
    {
        "qid": 6,
        "text1": "Find the pure liquid compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/j1JXocr.gif",
        "text2": "",
        "hints": ["The compound is soluble in water.", "Derivatives of this compound are used to create sturdy plastics."],
        "answers": ["propenoic acid", "acrylic acid", "CH2CHCOOH", "CH2=CHCOOH", "CH2CHCO2H", "CH2=CHCO2H"]
    },
    {
        "qid": 7,
        "text1": "Find the pure liquid compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/mi2jVRT.gif",
        "text2": "",
        "hints": ["Consider the molecular ion peak and nearby peaks. What do they tell you?", "The compound contains a halogen."],
        "answers": ["bromoacetic acid", "BrCH2COOH", "BrCH2CO2H", "CH2BrCOOH", "CH2BrCO2H"]
    },
    {
        "qid": 8,
        "text1": "Find the pure liquid compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/w1KE5uV.gif",
        "text2": "",
        "hints": ["Why are there no peaks between 30 and 127?", "The compound does not contain oxygen or nitrogen."],
        "answers": ["iodoethane", "CH3CH2I"]
    },
    {
        "qid": 9,
        "text1": "Find the gaseous compound that would give that mass spectrum.",
        "imgsrc": "https://i.imgur.com/Lr0jzrz.png",
        "text2": "",
        "hints": ["There are 2 constitutional isomers with the same functional group, both of which are valid answers to this question.", "The compound does not contain oxygen or nitrogen."],
        "answers": ["1-butene", "but-1-ene", "2-butene", "but-2-ene", "butene", "CH3CHCHCH3", "CH3CH=CHCH3", "CH2CH2CH2CH3", "CH2=CH2CH2CH3"]
    },
]

for question in questions:
    print(ms_qns.insert_one(question).inserted_id)