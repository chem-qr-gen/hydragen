import { now } from "jquery";

const AudioPlayer = {
  audio: new Audio(),
  isPlaying: false,
  load: function (url){
    this.audio = new Audio(url);
  },
  play: function () {
    this.audio.play();
    this.isPlaying = true;
  },
  pause: function () {
    this.audio.pause();
    this.isPlaying = false;
  },
  toggle: function () {
    if (this.audio.paused) {
      this.play();
    } else {
      this.pause();
    }
  },
  view: () => {}
};

const musicList = {
  1: {
    1: 'guitar_guitar.wav',
    2: 'guitar_oboe.wav',
    3: 'guitar_violin.wav',
  },
  2: {
    1: 'oboe_guitar.wav',
    2: 'oboe_oboe.wav',
    3: 'oboe_violin.wav',
  },
  3: {
    1: 'violin_guitar.wav',
    2: 'violin_oboe.wav',
    3: 'violin_violin.wav',
  },
}

function randomInt(min, max) {
  max++; // to make max inclusive
  const minCeiled = Math.ceil(min);
  const maxFloored = Math.floor(max);
  return Math.floor(Math.random() * (maxFloored - minCeiled) + minCeiled); // The maximum is exclusive and the minimum is inclusive
}

function randomMusic() {
  var nowMusic = randomInt(1, 3);
  var nextMusic = randomInt(1, 3);
  console.log("nowMusic: " + nowMusic + ", nextMusic: " + nextMusic);
  console.log(musicList);
  console.log(musicList[nowMusic][nextMusic]);
  var music = {
    playMusic: musicList[nowMusic][nextMusic],
    nextInstrument: nextMusic,
  };
  return music;
}

function nextMusic(nextInstrument) {
  var nowMusicNext = randomInt(1, 3);
  console.log(nowMusicNext);
  var music = {
    playMusic: musicList[nextInstrument][nowMusicNext],
    nextInstrument: nowMusicNext,
  };
  return music;
}

export {
  AudioPlayer,
  randomMusic,
  musicList,
  nextMusic
};

